/*
 * ReducedMesh.hpp
 *
 *  Created on: Dec 10, 2015
 *      Author: DAVIDE FORTI
 *        Mail: davide.forti@epfl.ch
 */

#ifndef REDUCEDMESH_HPP_
#define REDUCEDMESH_HPP_

namespace LifeV
{
//! ReducedMesh.hpp
/*!
    @author
    Reduced mesh

 */
template<typename MeshType>
class ReducedMesh
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                             mesh_Type;
    typedef boost::shared_ptr<MeshType>          meshPtr_Type;

    ReducedMesh( boost::shared_ptr< Epetra_Comm > comm, meshPtr_Type _mesh, int flag );

    ~ReducedMesh();

private:

    void countElementPerFlag();
    void allocatePerElements();
    void fillElementPerFlag();

public:

    void printElementPerFlag();
    void makeSubDivision();
    const UInt getNumElements() const {return M_numElementPerFlag; };
    const UInt *getSubmesh() const { return M_elements;};
    const UInt getFlag() const {return M_myflag; };

private:

    bool M_verbose;
    boost::shared_ptr< Epetra_Comm > M_comm;

    UInt *                           M_elements;
    int     						 M_numElementPerFlag;
    meshPtr_Type                     M_mesh;
    int 							 M_myflag;

};


template<typename MeshType>
ReducedMesh<MeshType>::ReducedMesh( boost::shared_ptr< Epetra_Comm > comm, meshPtr_Type _mesh, int flag ) :
	M_comm( comm ),
	M_mesh( _mesh ),
	M_myflag( flag ),
	M_numElementPerFlag ( 0 )
{}

template<typename MeshType>
ReducedMesh<MeshType>::~ReducedMesh()
{
    delete [] M_elements;
}

template<typename MeshType>
void ReducedMesh<MeshType>::countElementPerFlag()
{
    UInt nbElements ( M_mesh->numElements() );

    for (UInt iElement(0); iElement < nbElements; iElement++)
    {
        UInt markerID = M_mesh->element( iElement ).markerID( );

        if( M_myflag == markerID )
        {
            ++M_numElementPerFlag;
        }
    }
}

template<typename MeshType>
void ReducedMesh<MeshType>::fillElementPerFlag()
{
	int counter = 0;
    UInt nbElements( M_mesh->numElements( ) );

    for (UInt iElement(0); iElement < nbElements; iElement++)
    {
        // Extracting the marker
        UInt markerID = M_mesh->element( iElement ).markerID( );

        if( M_myflag == markerID )
        {
        	M_elements[ counter ] = iElement;
        	counter = counter + 1;
        }
    }
}

template<typename MeshType>
void ReducedMesh<MeshType>::allocatePerElements()
{
	M_elements = new UInt[M_numElementPerFlag];
}

template<typename MeshType>
void ReducedMesh<MeshType>::printElementPerFlag()
{
	std::cout << "ID " << M_comm->MyPID()
              << " flag " << M_myflag
              << " numElements: " << M_numElementPerFlag
              << std::endl;
}

template<typename MeshType>
void ReducedMesh<MeshType>::makeSubDivision()
{
    countElementPerFlag();
	allocatePerElements();
    fillElementPerFlag();

}

} // end LfeV namespace

#endif /* REDUCEDMESH_HPP_ */
