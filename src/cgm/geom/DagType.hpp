#ifndef DAG_TYPE_HPP
#define DAG_TYPE_HPP



class DagType
{

  public:
    
    enum FunctionalType {
      BasicTopologyEntity_TYPE = 0,
      SenseEntity_TYPE = 1,
      GroupingEntity_TYPE = 2 };
  
      // Constructor
      // 
      // Construct invalid DagType
    inline DagType()
      : intType(-1)
      {}
    
      // Constructor
      //
      // GroupingEntities and SenseEntities have the same 
      // dimension as their child BasicTopologyEntity type.
      // Dimension may be outside the range [0,3], but
      // if so is_valid() will return false.  (i.e. it is
      // possible to construct invalid DagType objects.)
    inline DagType( int my_dimension, FunctionalType func )
      : intType( 3 * my_dimension + func ) 
      { }
    
    
      // Return the dimension of this type.  
      //
      // GroupingEntities and SenseEntities have the same 
      // dimension as their child BasicTopologyEntity type.
    inline int dimension() const
      { return intType / 3; }
      
      // Return the base type of this subtype.
      // 
      // Return values are one of GroupingEntity_TYPE,
      // SenseEntity_TYPE, or BasicTopologyEntity_TYPE.
    inline FunctionalType functional_type() const
      { return static_cast<FunctionalType> (intType % 3); }
      
    
      // Is this a valid DagType?  
      //
      // An invalid type can be obtained by constructing
      // a DagType object with a dimension not in the 
      // range [0,3], or using the operators below to
      // move above a Body or below a RefVertex in the DAG.
    inline bool is_valid() const
      { 
        const unsigned my_body_type = 11;
        return static_cast<unsigned>(intType) <= my_body_type;
      }
      
      // get parent of this type
    inline DagType parent() const
      { return DagType(intType + 1); }
    
      
      // Comparison Operators.
      //
      //   An entity higher in the DAG (closer to Body) is
      //   considered greater than an entity lower in the
      //   DAG (closer to RefVertex).    

    inline bool operator==( DagType other ) const
      { return intType == other.intType; }

    inline bool operator!=( DagType other ) const
      { return intType != other.intType; }

    inline bool operator<=( DagType other ) const
      { return intType <= other.intType; }

    inline bool operator>=( DagType other ) const
      { return intType >= other.intType; }

    inline bool operator< ( DagType other ) const
      { return intType < other.intType; }

    inline bool operator> ( DagType other ) const
      { return intType > other.intType; }
    

      // Addition operators.
      //
      //  Increasing the value of the type moves the
      //  type up the DAG (closer to Body).  No bounds
      //  checking is done on the value.  Moving upwards
      //  from Body or downwards past RefVertex results
      //  in a DagType for which is_valid() returns false.

    inline DagType operator++()
      { return DagType(++intType); }

    inline DagType operator--()
      { return DagType(--intType); }

    inline DagType operator++(int)
      { return DagType(intType++); }

    inline DagType operator--(int)
      { return DagType(intType--); }

    inline DagType operator+=( int i )
      { return DagType(intType += i); }

    inline DagType operator-=( int i )
      { return DagType(intType -= i); }
      
    inline DagType operator+( int i ) const
      { return DagType(intType + i); }
      
    inline DagType operator-( int i ) const
      { return DagType(intType - i); }


      // Cast to bool.
      //
      // Same as is_valid().
    inline operator bool() const
      { return is_valid(); }
  

      // Get distance between types in DAG.
      //
      // Return value is negative if the passed type
      // is higher (closer to Body) in the DAG than
      // this type is.  The return value is 1 (one) 
      // if the passed type is an immediate child of 
      // this type.  The return value is 0 (zero) if
      // the types are the same.
    inline int operator-( DagType other ) const
      { return intType - other.intType; } 
  
      // Constants
inline static DagType       body_type() {return DagType(3,GroupingEntity_TYPE);}
inline static DagType  co_volume_type() {return DagType(3,SenseEntity_TYPE);}
inline static DagType ref_volume_type() {return DagType(3,BasicTopologyEntity_TYPE);} 
inline static DagType      shell_type() {return DagType(2,GroupingEntity_TYPE);}
inline static DagType    co_face_type() {return DagType(2,SenseEntity_TYPE);}
inline static DagType   ref_face_type() {return DagType(2,BasicTopologyEntity_TYPE);}
inline static DagType       loop_type() {return DagType(1,GroupingEntity_TYPE);}
inline static DagType    co_edge_type() {return DagType(1,SenseEntity_TYPE);}
inline static DagType   ref_edge_type() {return DagType(1,BasicTopologyEntity_TYPE);}
inline static DagType      chain_type() {return DagType(0,GroupingEntity_TYPE);}
inline static DagType  co_vertex_type() {return DagType(0,SenseEntity_TYPE);}
inline static DagType ref_vertex_type() {return DagType(0,BasicTopologyEntity_TYPE);}
inline static DagType    invalid_type() {return DagType(-1);}
  private:
  
    inline explicit DagType( int value )
      : intType( value ) {}
  
    int intType;
};


#endif

