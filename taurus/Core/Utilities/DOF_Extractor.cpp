#include <taurus/Core/Utilities/DOF_Extractor.hpp>
#include <chrono>
using namespace std::chrono;

namespace LifeV
{


DOF_Extractor::DOF_Extractor ( const fespacePtr_Type& fe_space ) :
	M_FESpace    ( fe_space ),
    M_status     ( 0 ),
    M_field      ( 1.0 ),
    M_dimension  ( 3.0 )
{
	M_markedDofs.reset( new VectorEpetra ( M_FESpace->map() ) );
}

void
DOF_Extractor::setup ( bcPtr_Type& bcHandler )
{
	M_markedDofs->zero();

	bcHandler->bcUpdate ( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );
	bcManageRhs ( *M_markedDofs, *M_FESpace->mesh(), M_FESpace->dof(), *bcHandler, M_FESpace->feBd(), 1.0, 1.0 );

	fullvector_indexes_unmarked.resize( M_FESpace->dof().numTotalDof()*3 );

	for ( int i = 0; i < M_markedDofs->epetraVector().MyLength(); ++i )
	{
		if ( (*M_markedDofs)[M_markedDofs->blockMap().GID(i)] == 1.0 )
		{
			indexes_marked.push_back ( M_markedDofs->blockMap().GID(i) );
		}
		else
		{
			fullvector_indexes_unmarked[M_markedDofs->blockMap().GID(i)] = 1.0;
			indexes_unmarked.push_back ( M_markedDofs->blockMap().GID(i) );
		}
	}

	fillMaps();
}

void
DOF_Extractor::fillMaps ( )
{
	buildMap ( indexes_marked, M_mapMarkedDofs_loc);

	buildMap ( indexes_unmarked, M_mapUnmarkedDofs_loc);
}


void
DOF_Extractor::buildMap ( std::vector<int> indices, mapPtr_Type& map)
{
	int* pointerToDofs (0);
	if (indices.size() > 0)
	{
		pointerToDofs = &indices[0];
	}

	map.reset ( new MapEpetra ( -1, static_cast<int> (indices.size() ), pointerToDofs, M_FESpace->map().commPtr() ) );

}

void
DOF_Extractor::restrictMatrixFast ( const matrixPtr_Type& matrix, const std::vector<int>& indexes, matrixPtr_Type& result)
{
	Real value = 1.0;
	matrixPtr_Type left;
	matrixPtr_Type right;
	std::vector<int>::const_iterator ITrow;

	//auto tStart = high_resolution_clock::now();
	left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 1 ) );
	left->zero();
	for (ITrow = indexes.begin(); ITrow != indexes.end(); ++ITrow)
	{
		if ( fullvector_indexes_unmarked[*ITrow] == 1.0 ) // test se e' interno
		{
			left->addToCoefficient( *ITrow , *ITrow, value );
		}
	}
	left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc );
	//auto tEnd = high_resolution_clock::now();
	//std::cout << "\n left build time = " << duration_cast<milliseconds>( tEnd - tStart ).count()/1000.0 << " s. ****\n";

	//right.reset ( new MatrixEpetra<Real> ( *left->transpose() ) );
	//right->globalAssemble(M_mapUnmarkedDofs_loc, M_FESpace->mapPtr() );

	//tStart = high_resolution_clock::now();
	right.reset ( new MatrixEpetra<Real> (*M_FESpace->mapPtr(), 1 ) );
	right->zero();
	for (ITrow = indexes.begin(); ITrow != indexes.end(); ++ITrow)
	{
		if ( fullvector_indexes_unmarked[*ITrow] == 1.0 ) // test se e' interno
		{
			right->addToCoefficient( *ITrow , *ITrow, value );
		}
	}

	right->globalAssemble(M_mapUnmarkedDofs_loc, M_FESpace->mapPtr() );
	//tEnd = high_resolution_clock::now();
	//	std::cout << "\n right build time = " << duration_cast<milliseconds>( tEnd - tStart ).count()/1000.0 << " s. ****\n";

	// DEVO FARE: result = left*matrix*(left').
	/*
	//Prima faccio tmp = matrix*(left'), poi result = left*tmp.
	matrixPtr_Type tmp ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );

	// tmp = matrix*(left')
	matrix->multiply(false, *right, false, *tmp, true);

	// result = left*tmp
	result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
	left->multiply(false, *tmp, false, *result, true);
	*/


	//Prima faccio tmp = left*matrix, poi result = tmp*right.
	matrixPtr_Type tmp ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );

	//tStart = high_resolution_clock::now();
	// tmp = left*matrix
	left->multiply(false, *matrix, false, *tmp, true);

	//tEnd = high_resolution_clock::now();
	//		std::cout << "\n left multiply time = " << duration_cast<milliseconds>( tEnd - tStart ).count()/1000.0 << " s. ****\n";

	//tStart = high_resolution_clock::now();

	// result = tmp*right
	result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
	tmp->multiply(false, *right, false, *result, true);

	//tEnd = high_resolution_clock::now();
	//	std::cout << "\n right multiply time = " << duration_cast<milliseconds>( tEnd - tStart ).count()/1000.0 << " s. ****\n";

}

void
DOF_Extractor::restrictMatrix ( const matrixPtr_Type& matrix, const std::string first, const std::string second, matrixPtr_Type& result)
{
	Real value = 1.0;
	matrixPtr_Type left;
	matrixPtr_Type right;
	std::vector<int>::const_iterator ITrow;
	std::vector<int>::const_iterator ITrow_col;

	if ( first.compare("unmarked") == 0 && second.compare("unmarked") == 0 )
	{
//		result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
//
////		for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
//		for (int z = 0; z < M_FESpace->map().map(Unique)->NumMyElements(); ++z)
//		{
//			if ( fullvector_indexes_unmarked[M_FESpace->map().map(Unique)->GID(z)] == 1.0 )
//			{
//				int numEntries;
//				int numEntriesPerRow;
//				int globalRow;
//
//				numEntries = matrix->matrixPtr()->NumMyEntries(z);
//				numEntriesPerRow = matrix->matrixPtr()->NumMyEntries(z);
//				globalRow = matrix->matrixPtr()->OperatorRangeMap().GID(z);
//
//				double* srcValues = new double[numEntries];
//				int* colIndices = new int[numEntries];
//
//				matrix->matrixPtr()->ExtractGlobalRowCopy(globalRow, numEntriesPerRow, numEntries, srcValues, colIndices);
//
//				std::vector<int> rows;
//				rows.push_back(globalRow);
//				std::vector<int> cols;
//				std::vector<Real> values;
//
//				for ( int i = 0; i < numEntriesPerRow; ++i )
//				{
//					//				if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
//					//				{
//					if ( fullvector_indexes_unmarked[colIndices[i]] == 1.0 )
//					{
//						cols.push_back(colIndices[i]);
//						values.push_back(srcValues[i]);
//					}
//					//				}
//				}
//
//				result->matrixPtr()->InsertGlobalValues( 1, &rows[0], cols.size(), &cols[0], &values[0] );
//
//				delete [] srcValues;
//				delete [] colIndices;
//			}
//		}
//		result->globalAssemble( M_mapUnmarkedDofs_loc, M_mapUnmarkedDofs_loc);
//		result->spy("matrice");


		left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 1 ) );

		for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
		{
//			if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
//			{
				left->addToCoefficient( *ITrow , *ITrow, value );
//				right->addToCoefficient( *ITrow , *ITrow, value );
//			}
		}

		left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc );
		right.reset ( new MatrixEpetra<Real> ( *left->transpose() ) );
		right->globalAssemble(M_mapUnmarkedDofs_loc, M_FESpace->mapPtr() );

		// DEVO FARE: result = left*matrix*(left'). Prima faccio tmp = matrix*(left'), poi result = left*tmp.

		matrixPtr_Type tmp ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );

		// tmp = matrix*(left')
		matrix->multiply(false, *right, false, *tmp, true);
		//tmp->globalAssemble( M_mapUnmarkedDofs_loc, M_FESpace->mapPtr() );

		// result = left*tmp
		result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
		//result->zero();

		left->multiply(false, *tmp, false, *result, true);
		//result->globalAssemble( M_mapUnmarkedDofs_loc, M_mapUnmarkedDofs_loc);
	}
	else if( first.compare("marked") == 0 && second.compare("marked") == 0 )
	{
		left.reset ( new MatrixEpetra<Real> (*M_mapMarkedDofs_loc, 50 ) );

		for (ITrow = indexes_marked.begin(); ITrow != indexes_marked.end(); ++ITrow)
		{
			if ( M_mapMarkedDofs_loc->map(Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
			{
				left->addToCoefficient( *ITrow , *ITrow, value );
			}
		}

		left->globalAssemble(M_FESpace->mapPtr(), M_mapMarkedDofs_loc);

		// DEVO FARE: result = left*matrix*(left'). Prima faccio tmp = matrix*(left'), poi result = left*tmp.

		matrixPtr_Type tmp ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );
		tmp->zero();

		// tmp = matrix*right
		matrix->multiply(false, *left, true, *tmp, false);
		tmp->globalAssemble( M_mapMarkedDofs_loc, M_FESpace->mapPtr() );

		// result = left*tmp
		result.reset ( new MatrixEpetra<Real> (*M_mapMarkedDofs_loc, 50 ) );
		result->zero();

		left->multiply(false, *tmp, false, *result, false);
		result->globalAssemble( M_mapMarkedDofs_loc, M_mapMarkedDofs_loc);
	}
	else if( first.compare("unmarked") == 0 && second.compare("marked") == 0 )
	{
		left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
		right.reset ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );

		for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
		{
			if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
			{
				left->addToCoefficient( *ITrow , *ITrow, value );
			}
		}

		for (ITrow = indexes_marked.begin(); ITrow != indexes_marked.end(); ++ITrow)
		{
			if ( M_mapMarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
			{
				right->addToCoefficient( *ITrow , *ITrow, value );
			}
		}


		left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc);
		right->globalAssemble(M_mapMarkedDofs_loc, M_FESpace->mapPtr());

		// DEVO FARE: result = left*matrix*(left'). Prima faccio tmp = matrix*(left'), poi result = left*tmp.

		matrixPtr_Type tmp ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );
		tmp->zero();

		// tmp = matrix*right
		matrix->multiply(false, *right, false, *tmp, false);
		tmp->globalAssemble( M_mapMarkedDofs_loc, M_FESpace->mapPtr() );

		// result = left*tmp
		result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
		result->zero();

		left->multiply(false, *tmp, false, *result, false);
		result->globalAssemble( M_mapMarkedDofs_loc, M_mapUnmarkedDofs_loc);
	}
}

void
DOF_Extractor::restrictVector ( const vectorPtr_Type& vector, const std::string part, vectorPtr_Type& result )
{
	Real value = 1.0;
	matrixPtr_Type left;
	std::vector<int>::const_iterator ITrow;

	if ( part.compare("marked") == 0 )
	{
		left.reset ( new MatrixEpetra<Real> (*M_mapMarkedDofs_loc, 50 ) );

		for (ITrow = indexes_marked.begin(); ITrow != indexes_marked.end(); ++ITrow)
		{
			if ( M_mapMarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
			{
				left->addToCoefficient( *ITrow , *ITrow, value );
			}
		}

		left->globalAssemble(M_FESpace->mapPtr(), M_mapMarkedDofs_loc);

		result. reset ( new VectorEpetra ( *M_mapMarkedDofs_loc, Unique ) );
		result->zero();

		*result += (*left * (*vector));
	}
	else if ( part.compare("unmarked") == 0 )
	{
		left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );

		for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
		{
			if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
			{
				left->addToCoefficient( *ITrow , *ITrow, value );
			}
		}


		left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc);

		result. reset ( new VectorEpetra ( *M_mapUnmarkedDofs_loc, Unique ) );
		result->zero();

		*result += (*left * (*vector));
	}
}

void
DOF_Extractor::restrictVectorFast ( const vectorPtr_Type& vector, const std::vector<int>& indexes, vectorPtr_Type& result )
{
	Real value = 1.0;
	matrixPtr_Type left;
	std::vector<int>::const_iterator ITrow;

	left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );

	for (ITrow = indexes.begin(); ITrow != indexes.end(); ++ITrow)
	{
		if ( fullvector_indexes_unmarked[*ITrow] == 1.0 ) // test se e' interno
		{
			left->addToCoefficient( *ITrow , *ITrow, value );
		}
	}

	left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc );

	result.reset ( new VectorEpetra ( *M_mapUnmarkedDofs_loc, Unique ) );
	result->zero();

	*result += (*left * (*vector));
}

void
DOF_Extractor::sumUnmarkedIntoGlobal( const vectorPtr_Type& input_vector, vectorPtr_Type& global_vector )
{
	Real value = 1.0;
	matrixPtr_Type left;
	std::vector<int>::const_iterator ITrow;

	left.reset ( new MatrixEpetra<Real> (M_FESpace->map(), 50 ) );

	for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
	{
		if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
		{
			left->addToCoefficient( *ITrow , *ITrow, value );
		}
	}

	left->globalAssemble(M_mapUnmarkedDofs_loc, M_FESpace->mapPtr());

	*global_vector += (*left * (*input_vector));
}

void
DOF_Extractor::restrictSnapshots ( const matrixPtr_Type& matrix, matrixPtr_Type& result )
{
	Real value = 1.0;
	matrixPtr_Type left;
	std::vector<int>::const_iterator ITrow;
	left.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );

	for (ITrow = indexes_unmarked.begin(); ITrow != indexes_unmarked.end(); ++ITrow)
	{
		if ( M_mapUnmarkedDofs_loc->map (Unique)->LID ( static_cast<EpetraInt_Type> (*ITrow) ) >= 0 )
		{
			left->addToCoefficient( *ITrow , *ITrow, value );
		}
	}

	left->globalAssemble(M_FESpace->mapPtr(), M_mapUnmarkedDofs_loc);

	result.reset ( new MatrixEpetra<Real> (*M_mapUnmarkedDofs_loc, 50 ) );
	result->zero();

	left->multiply(false, *matrix, false, *result, false);
	result->globalAssemble( matrix->domainMapPtr(), M_mapUnmarkedDofs_loc);
}

} // end namespace LifeV
