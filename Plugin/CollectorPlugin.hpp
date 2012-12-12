#ifndef COLLECTORPLUGIN_HPP
#define COLLECTORPLUGIN_HPP
#include "Vector3.hpp"


/*
 * CollectorPluginは、MolCollection -> Cell -> Mols -> SingleMolへと階層を下降しながら、指定された処理を実行するプラグイン。実行したい階層のDo*()をオーバーライドする。

 実行順序は保証されない。分子番号順に何かを実行させたい場合は別のPluginクラスを利用する。

 2010-2-9 例えば、外場を特定粒子だけに加えるのは、このプラグインでできるか?
 外場の処理はSystem/SimpleMD.cppのForce()の最後で行うべき。
 そこで活性化されるプラグインはない、と思う。
 対相互作用ならPairProcessPluginが使える。外場は対ではないから、このプラグインが近いと思うが、初期化(ファイルの読み込み)のあと、MDのコンストラクタを呼びだすところで、渡せるプラグインはProcessPluginだけじゃないかな。Force()から、ProcessPluginを呼びだせるようにするか。
 */

class MolCollection;
class Cell;
class Mols;
//class MonatomicMols;
class SingleMolEntity;


class CollectorPlugin{
public:
  virtual void HookL3( MolCollection* coll );
  virtual void HookL2( Cell* cell );
  virtual void HookL1( Mols* mols );
  virtual void HookL0( SingleMolEntity* mol ) = 0;
  //virtual void Hook( MonatomicMols* mols ) = 0;
};



class Test: public CollectorPlugin{
public:
  void HookL0( SingleMolEntity* mol );
};



//実用化テスト第一弾。
class QDoFPlugin: public CollectorPlugin{
  int dof;
public:
  QDoFPlugin(): dof(0) {}
  void HookL1( Mols* mol );
  void HookL0( SingleMolEntity* mol ){}
  int Result(){ return dof; }
};



//
class AddVelocityPlugin: public CollectorPlugin{
  Vector3 velocity;
public:
  AddVelocityPlugin( const Vector3& v ) : velocity(v) {}
  void HookL0( SingleMolEntity* mol );
  void HookL1( Mols* mols );
  virtual ~AddVelocityPlugin(){}
};



class ScaleVelocityPlugin: public CollectorPlugin{
  Vector3 ratio;
public:
  ScaleVelocityPlugin( const Vector3& r ) : ratio(r) {}
  void HookL0( SingleMolEntity* mol );
  virtual ~ScaleVelocityPlugin(){}
};



class ScaleVelocity2Plugin: public CollectorPlugin{
  double ratio;
public:
  ScaleVelocity2Plugin( double r ) : ratio(r) {}
  void HookL0( SingleMolEntity* mol );
  virtual ~ScaleVelocity2Plugin(){}
};



class ScalePositionPlugin: public CollectorPlugin{
  Vector3 ratio;
public:
  ScalePositionPlugin( const Vector3& r ) : ratio(r) {}
  void HookL0( SingleMolEntity* mol );
  virtual ~ScalePositionPlugin(){}
};




  

/*
//実用化テスト第二弾。なかなかうまくいかないもんだなあ。SerializeはCollectorの範疇ではない。
class SerializePlugin: public CollectorPlugin{
  double* p;
public:
  QDoFPlugin( double* p0 ): p(p0) {}
  void Hook( SingleMolEntity* mol );
};
*/


#endif
