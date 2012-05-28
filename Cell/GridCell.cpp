#include "Cell/GridCell.hpp"
#include "debug.hpp"


GridCell::GridCell()
{
  mesg("GridCell(): must be initialized explicitly.\n");
  exit(1);
}



GridCell::GridCell( const Box& box, int ncompo, int nx0, int ny0, int nz0 ) : boundary( box.clone() ), nx(nx0), ny(ny0), nz(nz0)
{
  int isP[3];
  int xmin = -1;
  int xmax = 1;
  int ymin = -1;
  int ymax = 1;
  int zmin = -1;
  int zmax = 1;
  if ( nx <= 3 ){
    nx = 1;
    xmin =xmax = 0;
  }
  if ( ny <= 3 ){
    ny = 1;
    ymin =ymax = 0;
  }
  if ( nz <= 3 ){
    nz = 1;
    zmin =zmax = 0;
  }
  nc  = nx*ny*nz;
  isP[0] = nx <= 3;
  isP[1] = ny <= 3;
  isP[2] = nz <= 3;
  bool isSemiperiodic = isP[0] || isP[1] || isP[2];

  //cells.resize(nc);
  Vector3 o;
  for( int iz=0; iz<nz; iz++ ){
    for( int iy=0; iy<ny; iy++ ){
      for( int ix=0; ix<nx; ix++ ){
	o.x = ix * box.x / nx;
	o.y = iy * box.y / ny;
	o.z = iz * box.z / nz;
	if ( isSemiperiodic ){
	  //����μ������Τ�PBC�ʥ���򰷤�����SemiperiodicBox�ˤ��롣
	  SemiperiodicBox smallbox(isP);
	  smallbox.Set( box.x / nx, box.y / ny, box.z / nz );
	  cells.push_back( RelocatableCellHandle( new RelocatableCell( smallbox, ncompo, o ) ) );
	}
	else{
	  Box smallbox;
	  smallbox.Set( box.x / nx, box.y / ny, box.z / nz );
	  cells.push_back( RelocatableCellHandle( new RelocatableCell( smallbox, ncompo, o ) ) );
	}
      }
    }
  }

  delta0 = new int [nc][13];
  //delta1 = new int [nc][26];
  //int c = 0;
  int d = 0;
  for( int z=zmin; z<=zmax; z++ ){
    int zz=z;
    if ( zz<0 ) zz += nz;
    for( int y=ymin; y<=ymax; y++ ){
      int yy=y;
      if ( yy<0 ) yy += ny;
      for( int x=xmin; x<=xmax; x++ ){
	int xx=x;
	if ( xx<0 ) xx += nx;
	//int p = (zz*ny+yy)*nx+xx;
	Vector3 o;

	//o���٤Υ��뤫�鸫�������Υ���ΰ��֡�������?
	o.x = -x * box.x / nx;
	o.y = -y * box.y / ny;
	o.z = -z * box.z / nz;
	if ( 0 < (z*ny+y)*nx+x ){
	  offset0[d] = o;
	  for( int iz=0; iz<nz; iz++ ){
	    int jz = iz+zz;
	    if ( nz <= jz ) jz -= nz;
	    for( int iy=0; iy<ny; iy++ ){
	      int jy = iy+yy;
	      if ( ny <= jy ) jy -= ny;
	      for( int ix=0; ix<nx; ix++ ){
		int jx = ix+xx;
		if ( nx <= jx ) jx -= nx;
		delta0[(iz*ny+iy)*nx+ix][d] = (jz*ny+jy)*nx+jx;
	      }
	    }
	  }
	  d++;
	  //cerr << x << " " << y << " " << z << endl;
	}
      }
    }
  }
  ndelta0 = d;
  //cerr << ndelta0 << " ndelta0\n";
}



GridCell::~GridCell()
{
  delete boundary;
  mesg("~GridCell\n");
  /*
    for( int i=0; i<nc; i++ ){
    delete cells->at(i);
    }
  */
  //free( cells );
}



void GridCell::SetBox( const Box& box )
{
  assert(0);
/*
  boundary = box;
  smallbox->x = box->x / nx;
  smallbox->y = box->y / ny;
  smallbox->z = box->z / nz;
  int c = 0;
  for( int z=0; z<nz; z++ ){
    for( int y=0; y<ny; y++ ){
      for( int x=0; x<nx; x++ ){
	Vector3 origin;
	origin.x   = x * smallbox->x;
	origin.y   = y * smallbox->y;
	origin.z   = z * smallbox->z;
	cells[c]->SetBox( smallbox, origin );
        c++;
      }
    }
  }
*/
}



void GridCell::Set( const int compo, const MolsHandle& mols )
{
  //First, initialize the subcells with the empty molshandle of the given molsentity
  for( int i=0; i<nc; i++ ){
    MolsHandle h = mols->EmptyClone( );
    cells[i]->Set( compo, h );
  }
  //distribute molecules to the cells
  int nmol = mols->Size();
  for( int i=0; i<nmol; i++ ){
    SingleMolHandle h = mols->Peek(i);
    h->BoxCoordinate( *boundary );
    Vector3 position = h->Position();
    //fprintf(stderr, "%f %f %f\n", position.x, position.y, position.z );
    int cell = whichcell( position );
    cells[cell]->PushAbs( compo, h );
  }
}



/*
void GridCell::ScaleVelocity( const Vector3& r )
{
  for( int i=0; i<nc; i++ ){
      cells[i]->ScaleVelocity( r );
  }
}
*/






void GridCell::ScaleRelativeCellVectors( const Vector3& r )
{
  //Scale relative cell vectors
  for( int i=0; i<ndelta0; i++ ){
    offset0[i].Scale( r );
  }
  //Ȣ����ɸ�����Х���٥��ȥ�򤽤줾��ۤʤ���ˡ�ǥ������뤷�Ƥ��롣r�Ϥ�Ȥ��1����ۤ�Τ鷺���������ʤ����ͤʤΤǡ����줾����ѿ��ؤθ������ѤΤ��������ۤʤ뤫�⤷��ʤ���
}





void GridCell::ProgressPosition( double dt )
{
  for( int i=0; i<nc; i++ ){
    cells[i]->ProgressPosition( dt );
  }
}



void GridCell::ProgressMomentum( double dt )
{
  for( int i=0; i<nc; i++ ){
    cells[i]->ProgressMomentum( dt );
  }
}



void GridCell::Preforce()
{
  //processes of order O(N)
  relocate();
  for( int i=0; i<nc; i++ ){
    cells[i]->Preforce();
  }
}



void GridCell::Postforce()
{
  //processes of order O(N)
  for( int i=0; i<nc; i++ ){
    if ( ! cells[i]->IsEmpty() )
      cells[i]->Postforce( );
  }
}



int GridCell::Force( const Combination& combi, const Truncation& rc, PotVir &pv )
{
  mesg("force0\n");
  //processes of order O(N^2)
  pv.Clear();
  int executed=0;
  //intra-cell
  for( int cell=0; cell<nc; cell++ ){
    PotVir p;
    //cerr << cell << " intra" << endl;
    if ( cells[cell]->Force( combi, rc, p ) ){
      executed ++;
      pv.Add(p);
    }
  }
  for( int cell=0; cell<nc; cell++ ){
    for( int nei=0; nei<ndelta0; nei++ ){
      PotVir p;
      int target = delta0[cell][nei];
      //cerr << cell << " " << nei << " " << target << " inter" << endl;
      //cerr << "(" << offset0[nei].x << "," << offset0[nei].y << "," << offset0[nei].z << ")" << endl;
      if ( cells[cell]->Force_Intercell( combi, offset0[nei], *(cells[target]), rc, p ) ){
        executed ++;
	pv.Add(p);
      }
    }
  }
  return executed;
}



int GridCell::PairProcess( PairProcessPlugin& p, const Combination& combi, const Truncation& rc ) const
{
  int executed=0;
  //intra-cell
  for( int cell=0; cell<nc; cell++ ){
    //cerr << cell << " intra" << endl;
    if ( cells[cell]->PairProcess( p, combi, rc ) ){
      executed ++;
    }
  }
  for( int cell=0; cell<nc; cell++ ){
    for( int nei=0; nei<ndelta0; nei++ ){
      int target = delta0[cell][nei];
      if ( cells[cell]->PairProcess_Intercell( p, combi, offset0[nei], *(cells[target]), rc ) ){
        executed ++;
      }
    }
  }
  return executed;
}



double GridCell::GetEk() const
{
  double ek=0;
  for( int i=0; i<nc; i++ ){
    ek += cells[i]->GetEk();
  }
  return ek;
}



void GridCell::GetEkt( Vector3& ekt ) const
{
  ekt.Set(0,0,0);
  for( int i=0; i<nc; i++ ){
    Vector3 ek;
    cells[i]->GetEkt( ek );
    ekt.x += ek.x;
    ekt.y += ek.y;
    ekt.z += ek.z;
  }
}



int GridCell::whichcell( const Vector3& coord )
{
  int x,y,z;
  //coord must be positive
  x = (int) (coord.x * nx / boundary->x );
  y = (int) (coord.y * ny / boundary->y );
  z = (int) (coord.z * nz / boundary->z );
  //fix error
  int err = 0;
  if ( nx <= x ){
    x -= nx;
    err += 1;
  }
  if ( ny <= y ){
    y -= ny;
    err += 2;
  }
  if ( nz <= z ){
    z -= nz;
    err += 4;
  }
  if ( err ){
    //printf("round %d\n", err );
  }
  return ((z*ny + y)*nx + x);
}



void GridCell::relocate()
{
  /*
   *MEMORY LEAK IS OCCURRING HERE.
   */
  for( int i=0; i<nc; i++ ){
    //get emmigrants
    RelocatableCell& from = *(cells[i]);
    MolsArray* em  = from.Emmigrate();
    MolsArray& emmigrants = *(em);
    //foreach component
    int ncompo = emmigrants.size();
    for( int compo=0; compo<ncompo; compo++ ){
      int nmol = emmigrants[compo]->Size();
      if (nmol){
	//foreach molecule
	for( int j=0; j<nmol; j++ ){
	  SingleMolHandle mol;
	  //peek the molecule
	  mol = emmigrants[compo]->Peek(j);
	  Vector3 position = mol->Position();
	  Vector3 offset(0,0,0);
          //cout << mol->GetOrder() << "\n";
	  //position must be in box coordinate
          //!!! it is not correct when the cell is semiperiodic.
          ///���Ȥ����������˥����ץ�ʾ��ˤϤ���Ϥޤ���!!!
          ///�����������ΤȤ��������������Ȥ�ͽ��Ϥʤ���
	  if ( position.x < 0 ) offset.x = boundary->x;
	  if ( position.y < 0 ) offset.y = boundary->y;
	  if ( position.z < 0 ) offset.z = boundary->z;
	  if ( boundary->x <= position.x ) offset.x = -boundary->x;
	  if ( boundary->y <= position.y ) offset.y = -boundary->y;
	  if ( boundary->z <= position.z ) offset.z = -boundary->z;
	  mol->Translate( offset );
	  position.x += offset.x;
	  position.y += offset.y;
	  position.z += offset.z;
	  //assing the cell
	  int cell = whichcell( position );
	  if ( i == cell ) {
	    cerr << "migration in the same cell." << mol->GetOrder()
		 << " " << position.x
		 << " " << position.y
		 << " " << position.z
		 << "\n";
	  }
	  //push into it
	  cells[cell]->PushAbs( compo, mol );
	}
      }
    }
    delete em;
  }
}



int GridCell::GetSize( const int compo ) const
{
  int size = 0;
  for( int i=0; i<nc; i++ ){
    size += cells[i]->GetSize( compo );
  }
  return size;
}



double
GridCell::GetMass() const
{
  double mass = 0;
  for( int i=0; i<nc; i++ ){
    mass += cells[i]->GetMass();
  }
  return mass;
}






int GridCell::GetNumCompo() const
{
  return cells[0]->GetNumCompo();
}




void GridCell::Write( const Unit& unit, ostream& to )
{
  to << "@VOXN\n"
     << nx
     << " " << ny
     << " " << nz << '\n';
  boundary->WriteBOX3( unit, to );
}




MolsHandle
GridCell::PullAll( const int compo )
{
  MolsHandle mh = cells[0]->mols[compo]->EmptyClone();
  //���ꤢ�碌�δؿ������Ǻ�äƤߤ롣
  for( int c=0; c<nc; c++ ){
    mh->Concat( cells[c]->PullAll( compo ) );
  }
  return mh;
}



MolsHandle
GridCell::CopyAll( const int compo )
{
  MolsHandle mh = cells[0]->mols[compo]->EmptyClone();
  //���ꤢ�碌�δؿ������Ǻ�äƤߤ롣
  for( int c=0; c<nc; c++ ){
    MolsHandle tmp = cells[c]->CopyAll( compo );
    tmp->Translate( cells[c]->GetOrigin() );
    mh->Concat( tmp );
  }
  return mh;
}



void GridCell::TotalMomentum( Vector3& momentum ) const
{
  momentum.Set(0,0,0);
  for( int c=0; c<nc; c++ ){
    Vector3 m;
    cells[c]->TotalMomentum( m );
    momentum.x += m.x;
    momentum.y += m.y;
    momentum.z += m.z;
  }
}



/*
void GridCell::AddVelocity( const Vector3& velocity )
{
  for( int c=0; c<nc; c++ ){
    cells[c]->AddVelocity( velocity );
  }
}
*/



void GridCell::Execute( CollectorPlugin& plugin )
{
  for( int c=0; c<nc; c++ ){
    cells[c]->Execute( plugin );
  }
}



void GridCell::Rescale( const Vector3& r )
{
  // GridCell�γ��Ȥ򥹥����뤹�롣
  boundary->Scale( r );
  // ��ξ�Ȣ�򥹥����뤹�롣
  for( int i=0; i<nc; i++ ){
      cells[i]->Rescale( r );
  }
  //���а��֥٥��ȥ�򥹥����뤹�롣
  ScaleRelativeCellVectors( r );
  //Ȣ����ɸ�����Х���٥��ȥ�򤽤줾��ۤʤ���ˡ�ǥ������뤷�Ƥ��롣r�Ϥ�Ȥ��1����ۤ�Τ鷺���������ʤ����ͤʤΤǡ����줾����ѿ��ؤθ������ѤΤ��������ۤʤ뤫�⤷��ʤ���
}



//��ɸ����Ф����ݥ��󥿤κ�������˽��֤�����롣
double*
GridCell::serialize( double* p ) const
{
  //�礭����������Ū�˽���������ö������serialize���Ƥ�餦��
  //��ʬ���ȤˤޤȤ��ɬ�פ�����Τ�����GridCell������ʬ�ϸ����ʤ��ΤǤ�?
  assert(0);
  /*
  double tmp[];
    for( int compo=0; compo < ncompo; compo++ ){
    for( int i=0; i<nc; i++ ){
      //cells[i]->Rescale( r );
    }
  }
  */
  return p;
}
