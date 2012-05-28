#ifndef EXCHANGEPT_HPP
#define EXCHANGEPT_HPP
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include "Random.hpp"
#include "Exchange.hpp"

using namespace boost;



class ExchangePT : public Exchange
{
public:
  ExchangePT( const RandomD01Handle& r );
  bool Trial( ExchangePT& target );
};


#endif
