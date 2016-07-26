#include <taurus/Core/ReducedBasis/ImportSnapshots.hpp>

namespace LifeV
{


ImportSnapshots::ImportSnapshots ( const fespacePtr_Type& fe_space, const int numParameters, const GetPot& dataFile, bcPtr_Type bc ) :
	M_FESpace        ( fe_space ),
	M_numParameters  ( numParameters ),
	M_datafile       ( dataFile )
{
	M_bc_markingDOFs = bc;
}

void
ImportSnapshots::readSnapshots( std::string field_name)
{
	buildColumnMap();

	M_solution.reset ( new VectorEpetra ( M_FESpace->map(), Unique ) );

	M_X.reset ( new MatrixEpetra<Real> ( M_FESpace->map(), M_numParameters) );

	std::string const fileName_importer =  M_datafile ( "importer/name_input_file", "laplaciano_ridotto");
	M_importer.reset ( new  ExporterHDF5<mesh_Type > ( M_datafile, fileName_importer) );
	M_importer->setMeshProcId ( M_FESpace->mesh(), M_FESpace->map().commPtr()->MyPID() );
	std::string iterationString;

	for (int i = 0; i < M_numParameters; ++i )
	{
		std::ostringstream iter;
		iter.fill ( '0' );
		iter << std::setw (5) << ( i );
		iterationString = iter.str();

		M_solution->zero();
		LifeV::ExporterData<mesh_Type> solutionReader (LifeV::ExporterData<mesh_Type>::ScalarField,
													   std::string (field_name + "." + iterationString),
													   M_FESpace,
													   M_solution,
													   UInt (0),
													   LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

		M_importer->readVariable ( solutionReader );

		for ( int j = 0; j < M_solution->epetraVector().MyLength(); ++j )
		{
			M_X->addToCoefficient( M_solution->blockMap().GID(j), i, (*M_solution)(M_solution->blockMap().GID(j)) );
		}
	}

	M_X->globalAssemble(M_map_column, M_FESpace->mapPtr());

	M_importer->closeFile();
}

void
ImportSnapshots::readSnapshotsVectorial( std::string field_name )
{
	buildColumnMap();

	M_solution.reset ( new VectorEpetra ( M_FESpace->map(), Unique ) );

	M_X.reset ( new MatrixEpetra<Real> ( M_FESpace->map(), M_numParameters) );

	std::string const fileName_importer =  M_datafile ( "importer/name_input_file", "laplaciano_ridotto");
	M_importer.reset ( new  ExporterHDF5<mesh_Type > ( M_datafile, fileName_importer) );
	M_importer->setMeshProcId ( M_FESpace->mesh(), M_FESpace->map().commPtr()->MyPID() );
	std::string iterationString;

	for (int i = 0; i < M_numParameters; ++i )
	{
		std::ostringstream iter;
		iter.fill ( '0' );
		iter << std::setw (5) << ( i );
		iterationString = iter.str();

		M_solution->zero();
		LifeV::ExporterData<mesh_Type> solutionReader ( LifeV::ExporterData<mesh_Type>::VectorField,
														std::string (field_name + "." + iterationString),
														M_FESpace,
														M_solution,
														UInt (0),
														LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

		M_importer->readVariable ( solutionReader );

		for ( int j = 0; j < M_solution->epetraVector().MyLength(); ++j )
		{
			M_X->addToCoefficient( M_solution->blockMap().GID(j), i, (*M_solution)(M_solution->blockMap().GID(j)) );
		}
	}

	M_X->globalAssemble(M_map_column, M_FESpace->mapPtr());

	M_importer->closeFile();
}

void
ImportSnapshots::buildColumnMap( )
{
	std::vector<int> column_indices;
	for ( int i = 0; i < M_numParameters; ++i )
		column_indices.push_back(i);

	int* pointerToDofs (0);
	if (column_indices.size() > 0)
		pointerToDofs = &column_indices[0];

	M_map_column.reset ( new MapEpetra ( -1, static_cast<int> (column_indices.size() ), pointerToDofs, M_FESpace->map().commPtr() ) );
}

void
ImportSnapshots::getMatrix ( matrixPtr_Type& matrix )
{
	M_bc_markingDOFs->bcUpdate( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );

	M_dofExtractor.reset ( new DOF_Extractor ( M_FESpace ) );

	M_dofExtractor->setup( M_bc_markingDOFs );

	M_dofExtractor->restrictSnapshots( M_X, matrix);
}

} // end namespace LifeV
