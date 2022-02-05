#ifndef CTRIPLE_H
#define CTRIPLE_H

#include <iostream>

//class Ctriple;
//bool Ctriple::operator<(Ctriple &a, Ctriple &b)

template <typename T1, typename T2, typename T3>
class Ctriple
{
 private:
  T1 t1;
  T2 t2;
  T3 t3;

 public:
 Ctriple(const T1 &pt1,const T2 &pt2,const T3 &pt3):t1(pt1),t2(pt2),t3(pt3)
  {}

  T1 first()  { return t1; }
  T2 second() { return t2; }
  T3 third()  { return t3; }

  friend bool operator<(const Ctriple<T1, T2, T3> &a, const Ctriple<T1, T2, T3> &b)
  {
    if (a.t1 != b.t1)
      return a.t1 < b.t1;
    if (a.t2 != b.t2)
      return a.t2 < b.t2;
    return a.t3 < b.t3;
  }

  void print(std::ostream &os) const
  {
    os << "(" << t1 << "," << t2 << "," << t3 << ")";
  }

};


template <typename T1, typename T2, typename T3>
inline std::ostream &operator<<(std::ostream &os, const Ctriple<T1, T2, T3> &t)
{
  t.print(os);
  return os;
}









#endif


