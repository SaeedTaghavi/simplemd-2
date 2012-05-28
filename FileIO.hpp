#ifndef FILEIO_HPP
#define FILEIO_HPP

#include <cstdio>
#include <map>
#include "Unit.hpp"
#include "MolProperty.hpp"
#include "SingleMol/RigidBody2.hpp"
#include "Mols/RigidBodies2.hpp"
#include "Mols/RigidBodies.hpp"

FlexibleHandle
ReadDEFR( string id08, const Unit& u, FILE* file );
FlexibleHandle
ReadDEFP( string id08, const Unit& u, FILE* file );
void
WriteMDVW( const RigidBody2& r, const Unit& unit, ostream& to );
void
WriteMDVWHeader( const MolPropertyHandle& mph, const RigidBodies2& r, const BoxHandle& box, const Unit& unit, ostream& to );
void
WriteMDVW( RigidBodies2& r, const Unit& unit, ostream& to );
void
Read( PolyatomicMol& m, int o, FILE* file, const Unit& u );
void
WriteDTPS( double dt, const Unit& unit, ostream& to );
void
Write( map<string, FlexibleHandle>& moldict, const Unit& unit, ostream& to );
double* serializeforce0( const MolsHandle m, double* ptr );
double* unserialize0( MolsHandle m, double* ptr );


#endif
