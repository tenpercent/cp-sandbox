///---------------------------------------------------------------------
/// Class: TDVector.hpp
/// Description: Holds some data for the DecompSweepTool.
/// Owner: David R. White
/// Creation Date: 2/25/2003
///---------------------------------------------------------------------

#ifndef TDVECTOR_HPP
#define TDVECTOR_HPP

#include "ToolData.hpp"
#include "CubitVector.hpp"
#include "CastTo.hpp"

class TDVector : public ToolData
{
private:
  
  const CubitVector originalPosition;
  
public:
  
  TDVector(const CubitVector &position)
    : originalPosition(position) {}
    ///
    ///constructor
    ///

  ~TDVector()
    {}
    ///
    /// Destructor
    ///

  const CubitVector& get_vector() const
  {
    return originalPosition;
  }
    ///
    /// Return the original position.
    ///
      
  static int is_td_vector(const ToolData* td)
  {
    return (CAST_TO(td, const TDVector) != NULL );
  }

};


#endif // TDVector

