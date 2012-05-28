#include "Cell/PeriodicCell.hpp"

PeriodicCell::PeriodicCell( const BoxHandle& box, int ncompo ) : SimpleCell( ncompo )
{
  boundary = box;
}



PeriodicCell::~PeriodicCell()
{
  mesg("~PeriodicCell\n");
}



void PeriodicCell::Preforce()
{
  relocate();
  SimpleCell::Preforce();
}



int PeriodicCell::Force( const Combination& combi, const Truncation& rc, PotVir &pv )
{
  pv.ep = 0;
  pv.vr.Clear();
  int executed = 0;
  int ncombi = combi.ima.size();
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    PotVir p;
    InteractionMatrix* im = combi.ima[co].get();
    int result;
    if ( compo1 == compo2 ){
      result = mols[compo1]->Force_PBC( *im, boundary, rc, p );
    }
    else{
      result = mols[compo1]->Force_PBC( mols[compo2], *im, boundary, rc, p );
    }
    if ( result ){
      executed ++;
      printf("%d %d %f\n", compo1, compo2, p.ep );
      pv.ep += p.ep;
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          pv.vr.v[i][j] += p.vr.v[i][j];
        }
      }
    }
  }
  return executed;
}



void PeriodicCell::relocate()
{
  int ncompo = mols.size();
  for( int i=0;i<ncompo;i++ ){
    mols[i]->BoxCoordinate( boundary );
  }
}



BookPeriodicCell::BookPeriodicCell( const BoxHandle& box, int ncompo ) : PeriodicCell( box, ncompo )
{
  //preset simplepair
  simplepair = vector<SimplePair>( ncompo*(ncompo+1)/2 );
}



int BookPeriodicCell::Force( const Combination& combi, const Truncation& rc, PotVir &pv )
{
  pv.ep = 0;
  pv.vr.Clear();
  int executed = 0;
  int ncombi = combi.ima.size();
  //
  //���åȥ��յ�Υ����1AĹ��ǲ��Хꥹ�Ȥ�������롣
  //
  for( int co=0; co<ncombi; co++ ){
    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    const Mols& m1 = *(mols[compo1]);
    const Mols& m2 = *(mols[compo2]);
    //
    //�����ϡ����Ū��book��ƹ��ۤ��ʤ���Ф����ʤ�����
    //���ΤȤ������Τ�book���ۤ�Ԥ����������äơ�Ĺ���֤η׻��ϤǤ��ʤ���
    //benchmark�ѡ�
    //
    if ( simplepair[co].listitems.size() == 0 ){
      int npair;
      if ( compo1 == compo2 ){
	npair = simplepair[co].list_pbc( m1, boundary, rc.outer+1.0 );
      }
      else{
	npair = simplepair[co].list_pbc( m1, m2, boundary, rc.outer+1.0 );
      }
      fprintf(stderr, "Raw pairs(%d, %d) = %d\n", compo1, compo2, npair);
    }
  }
  //
  //���Хꥹ�Ȥΰ��̡��ºݤ���ߺ��Ѥ����ФΤߤ���Ф��롣������Ϸ׻��˲󤹡�
  //
  for( int co=0; co<ncombi; co++ ){

    int compo1 = combi.compo1[co];
    int compo2 = combi.compo2[co];
    const Mols& m1 = *(mols[compo1]);
    const Mols& m2 = *(mols[compo2]);
    PotVir p;
    InteractionMatrix* im = combi.ima[co].get();
    int result;
    TruncPair* truncpair;
    truncpair = simplepair[co].Compress( m1, m2, boundary, rc );
    fprintf(stderr, "Compressed pairs(%d, %d) = %d\n", compo1, compo2, truncpair->listitems.size());
    result = mols[compo1]->Force_General( mols[compo2], *im, *truncpair, p );
    delete truncpair;
    if ( result ){
      executed ++;
      printf("%d %d %f\n", compo1, compo2, p.ep );
      pv.ep += p.ep;
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          pv.vr.v[i][j] += p.vr.v[i][j];
        }
      }
    }
  }
  return executed;
}



