#ifndef SYSTEM_HPP
#define SYSTEM_HPP

//#include "Vector3.hpp"

class Vector3;

/*
 *Integratorの種類にあわせ、MDの進め方が違ってくる。その違いを吸収する
 *ためのクラス。
 *
 *MDのメインプログラムからは、ここのクラスのひとつを選択し、
 *初期化フェーズではinitialize()を、
 *実行フェーズではrunnungを、
 *終了フェーズではterminateを呼びだす。
 *
 *MDの過程で、他の処理をやらせたければ、これらのクラスを継承するか、包
 *含した新クラスを作ればよい。
 */



class Main
{
public:
  virtual void initialize() = 0;
  virtual bool running() = 0;
  virtual void terminate() = 0;
};






class VirtualMD
{
public:
  virtual void initialize() = 0;
  virtual bool running() = 0;
  virtual void terminate() = 0;
  virtual double OneStep() = 0;
  virtual double GetEk() = 0;
  virtual void GetEkt( Vector3& ) = 0;
};





#endif
