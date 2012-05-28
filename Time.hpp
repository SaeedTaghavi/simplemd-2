#ifndef TIME_HPP
#define TIME_HPP

/*
 *timekeeper module
 */

class TimeKeeper{
private:
  double absolute;
  double delta;
public:
  TimeKeeper( double abs, double del )
  {
    Absolute( abs );
    Delta( del );
  }
  void Absolute( double abs )
  {
    absolute = abs;
  }
  void Delta( double del )
  {
    delta = del;
  }
  double Absolute() const
  {
    return absolute;
  }
  double Delta() const
  {
    return delta;
  }
  double Progress( double ratio )
  {
    absolute += ratio * delta;
    return absolute;
  }
};





#endif
