

#include "DagDrawingTool.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "RefVolume.hpp"
#include "Body.hpp"
#include "DLIList.hpp"
#include "CubitColorConstants.hpp"
#include "SenseEntity.hpp"
#include "CastTo.hpp"
#include "CoVertex.hpp"
#include "Chain.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "CoFace.hpp"
#include "Shell.hpp"

//constants:


const char* const DDT_BODY_ABBREV    = "B";
const char* const DDT_COVOL_ABBREV   = "CVm";
const char* const DDT_REFVOL_ABBREV  = "Vm";
const char* const DDT_SHELL_ABBREV   = "Sh";
const char* const DDT_COFACE_ABBREV  = "CF";
const char* const DDT_REFFACE_ABBREV = "F";
const char* const DDT_LOOP_ABBREV    = "L";
const char* const DDT_COEDGE_ABBREV  = "CE";
const char* const DDT_REFEDGE_ABBREV = "E";
const char* const DDT_CHAIN_ABBREV   = "Ch";
const char* const DDT_COVTX_ABBREV   = "CVx";
const char* const DDT_REFVTX_ABBREV  = "Vx";
//Abbreviations used for type names.
//Note: If abbreviations are not used, the
// full type name is taken from 
// CubitEntity::entity_name( EntityType ).


//const CubitLink DDT_DEFAULT_TRAVERSE_TYPE = CUBIT_ALL_LINKS;

DagDrawingTool* DagDrawingTool::instance_ = 0;

//constats for query direction
const int DDT_UP_DAG   =  1;
const int DDT_DOWN_DAG = -1;

//bit masks for link type
const int DDT_LINK_NONE = 0;
const int DDT_LINK_UP   = 1;
const int DDT_LINK_DOWN = 2;
const int DDT_LINK_BOTH = 3;
const int DDT_BAD_LINK  = 4;

DagDrawingTool::DagDrawingTool()
{
	node_table = 0;
}

void DagDrawingTool::one_sort_pass( int start, int end, float* val_list )
{
	int i, j, k;
	DagNodeTable& table = *node_table;
	//int row_count = table.rows();
	int cur_row_len = 0;
	if( start == end ) return;
	int dir = (start < end) ? 1 : -1;
	
	//for each parent row of the passed row of nodes...
	for( i = start; i != (end + dir); i += dir )
	{
		//populate the array with the values to sort by
		cur_row_len = table[i].length();
		for( j = 0; j < cur_row_len; j++ )
		{
			val_list[j] = 0.0;
			int link_count = 0;
			for( k = 0; k < table[i-dir].length(); k++ )
			{
				if( link_type( i, j, i - dir, k ) )
				{
					val_list[j] += (float)k;
					link_count++;
				}
			}
			if( link_count ) val_list[j] /= (float)link_count;
		}		
		
		//sort the row by the values in value_list
		for( j = 0; j < cur_row_len - 1; j++ )
		{
			int smallest = j;
			for( k = j + 1; k < cur_row_len; k++ )
			{
				if( val_list[k] < val_list[smallest] )
					smallest = k;
			}
			if( smallest != j )
			{
                                  //float temp = val_list[j];
				val_list[j] = val_list[smallest];
				val_list[smallest] = val_list[j];
				
				ModelEntity* node = table[i][j];
				table[i][j] = table[i][smallest];
				table[i][smallest] = node;
			}
		}
	}
	
}	
/*
void DagDrawingTool::sort_table( )
{
	DagNodeTable& table = *node_table;
	int passes = table.rows() - 1;
	int i;
	
	int max_row_size = 0;
	int row_count = table.rows();
	for( i = 0; i < row_count; i++ )
		if( table[i].length() > max_row_size )
			max_row_size = table[i].length();
	float* value_list = new float[max_row_size];

	for( i = 0; i < passes; i++ )
	{
		one_sort_pass( 1, row_count - 1, value_list );
		one_sort_pass( row_count - 2, 0, value_list );
	}
		
	delete [] value_list;
}
*/ 

void DagDrawingTool::sort_row( int row, int wrt_dir )
{
  assert((wrt_dir == 1) || (wrt_dir == -1));
  
  DagNodeTable& table = *node_table;
  int row_len = table[row].length();
  int rel_row = row + wrt_dir;
  if( (row_len < 2) || (rel_row < 0) || (rel_row == table.rows()) )
    return;
  
  int* index_array = new int[row_len];
  int num_rels = table[rel_row].length();
  int i, j, k, l, m;
  
  //Append first node to list
  index_array[0] = 0;
  
  //For each remaining node, insert it in the location at which
  //it causes the least intersections of links.
  for( i = 1; i < row_len; i++ )
  {
    int best_index = 0;
    int intersection_count = 0;
    
    //Get the intersection count if the node were inserted first
    //in the row.

    //For each possible parent node
    for( k = 1; k < num_rels; k++ )
    {
      //Is it a parent?
      if( ! link_type( row, i, rel_row, k ) )
        continue;
      
      //Count how many other links that link intersects.
      
      //For all the other nodes in the list...
      for( l = 0; l < i; l++ )
      {
        //For each parent to the left of my parent
        for( m = 0; m < k; m++ )
        {
          if( link_type( row, index_array[l], rel_row, m ) )
            intersection_count++;
        }
      }
    }
    
    //Now check the intersection count for each additional
    //possible location of insertion.       
    for( j = 1; j <= i; j++ )
    {
      int current_count = 0;
      
      //For each possible parent node
      for( k = 0; k < num_rels; k++ )
      {
        //Is it a parent?
        if( ! link_type( row, i, rel_row, k ) )
          continue;
          
        //Count how many other links from nodes to the
        //left of me intersect with my link to this parent.
        for( l = 0; l < j; l++ )  
        {
          //For each parent to the right of my parent
          for( m = k+1; m < num_rels; m++ )
          {
            if( link_type( row, index_array[l], rel_row, m ) )
              current_count++;
          }
        }
        
        //Count how many other links from nodes to the
        //right of me intersect with my link to this parent.
        for( l = j; l < i; l++ )
        {
          //For each parent to the left of my parent
          for( m = 0; m < k; m++ )
          {
            if( link_type( row, index_array[l], rel_row, m ) )
              current_count++;
          }
        }
      }
      
      if( current_count <= intersection_count )
      {
        intersection_count = current_count;
        best_index = j;
      }
    }
    
    //Now insert this node at the best location.
    for( j = i; j > best_index; j-- )
      index_array[j] = index_array[j-1];
    index_array[best_index] = i;
  }
  
  //Now rearrange the nodes in the actual row to match 
  //what is stored in index_list
  ModelEntity** node_array = new ModelEntity*[row_len];
  for( i = 0; i < row_len; i++ )
  {
    node_array[i] = table[row][index_array[i]];
  }
  for( i = 0; i < row_len; i++ )
  {
    table[row][i] = node_array[i];
  }
  delete [] node_array;
  delete [] index_array;
}
        
        
/*
void DagDrawingTool::draw_table()
{
	int i, j, k;
	DagNodeTable& table = *node_table;
	for( i = 0; i < table.rows(); i++ )
	{
		int len = table[i].length();
		for( j = 0; j < len; j++ )
		{
			draw_node( i, j );
			if( i > 0 ) 
			{
				int len2 = table[i-1].length();
				for( k = 0; k < len2; k++ ) draw_link( i, j, i - 1, k );
			}
      if( draw_host_para_links_ )
      {
        for( k = 0; k < j; k++ )
          draw_host_link( i, k, j );
		  }
    }
	}
}


void DagDrawingTool::draw_DAG( DLIList<ModelEntity*>& node_list, int up, int down )
{
	int i;
	assert( (up >= 0) && (down >= 0) );
	if( node_list.size() == 0 ) return;
	
	DLIList<ModelEntity*> relatives[2];
	int num_rows = up + down + 1;
	
	assert( ! node_table );
	node_table = new DagNodeTable( num_rows );
	start_row = down;
	node_table->initialize_row( start_row, node_list );
	
	relatives[ start_row % 2 ] = node_list;
	for( i = start_row - 1; i >= 0; i-- )
	{
		relatives[i % 2].clean_out();
		get_relatives( relatives[(i + 1) % 2], relatives[i % 2], DDT_DOWN_DAG );
		node_table->initialize_row( i, relatives[i % 2] );
	}
	relatives[ start_row % 2 ] = node_list;
	for( i = start_row + 1; i < num_rows; i++ )
	{
		relatives[i % 2].clean_out();
		get_relatives( relatives[(i - 1) % 2], relatives[i % 2], DDT_UP_DAG );
		node_table->initialize_row( i, relatives[i % 2] );
	}
	

	sort_table();

	initialize_graphics();
	draw_table();
	finalize_graphics();

	delete node_table;
	node_table = 0;
}
*/


DagDrawingTool* DagDrawingTool::instance()
{
	return instance_ ? instance_ : instance_ = new DagDrawingTool();
}

DagDrawingTool::~DagDrawingTool()
{
	instance_ = 0;
	assert( ! node_table );
}

/*
void DagDrawingTool::draw_DAG( Body*      body,           int down ) 
{
  DLIList<RefEntity*> entity_list;
  entity_list.append(body);
  draw_DAG( entity_list, 0, down );
}

void DagDrawingTool::draw_DAG( RefVolume* volume, int up, int down )
{
  DLIList<RefEntity*> entity_list;
  entity_list.append(volume);
  draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( RefFace*   face,   int up, int down )
{
  DLIList<RefEntity*> entity_list;
  entity_list.append(face);
  draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( RefEdge*   edge,   int up, int down )
{
  DLIList<RefEntity*> entity_list;
  entity_list.append(edge);
  draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( RefVertex* vertex, int up           )
{
  DLIList<RefEntity*> entity_list;
  entity_list.append(vertex);
  draw_DAG( entity_list, up, 0 );
}

void DagDrawingTool::draw_DAG( DLIList<RefEdge*>& edge_list, int up, int down )
{
	DLIList<RefEntity*> entity_list;
	CAST_LIST_TO_PARENT( edge_list, entity_list );
	draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( DLIList<RefFace*>& face_list, int up, int down )
{
	DLIList<RefEntity*> entity_list;
	CAST_LIST_TO_PARENT( face_list, entity_list );
	draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( DLIList<RefVertex*>& vtx_list, int up )
{
	DLIList<RefEntity*> entity_list;
	CAST_LIST_TO_PARENT( vtx_list, entity_list );
	draw_DAG( entity_list, up, 0 );
}

void DagDrawingTool::draw_DAG( DLIList<RefVolume*>& vol_list, int up, int down )
{
	DLIList<RefEntity*> entity_list;
	CAST_LIST_TO_PARENT( vol_list, entity_list );
	draw_DAG( entity_list, up, down );
}

void DagDrawingTool::draw_DAG( DLIList<Body*>& body_list, int down )
{
	DLIList<RefEntity*> entity_list;
	CAST_LIST_TO_PARENT( body_list, entity_list );
	draw_DAG( entity_list, 0, down );
}

void DagDrawingTool::draw_DAG( DLIList<RefEntity*>& entity_list, int up, int down )
{
	DLIList<ModelEntity*> node_list;
	for( int i = entity_list.size(); i > 0; i-- )
	{
		RefEntity* re_ptr = entity_list.get_and_step();
		ModelEntity* me_ptr = CAST_TO( re_ptr, ModelEntity);
		node_list.append_unique( me_ptr );
	}
	draw_DAG( node_list, up, down );
}
*/


/*
void DagDrawingTool::draw_host_parasite_links( CubitBoolean yn )
{ draw_host_para_links_ = yn; }
CubitBoolean DagDrawingTool::draw_host_parasite_links() const
{ return draw_host_para_links_; }
void DagDrawingTool::primary_link_color( int color )
{ prim_link_color_ = color; }
void DagDrawingTool::secondary_link_color( int color )
{ sec_link_color_ = color; }
void DagDrawingTool::host_para_link_color( int color )
{ host_para_color_ = color; }
void DagDrawingTool::bad_link_color( int color )
{ bad_link_color_ = color; }
void DagDrawingTool::node_color( int color )
{ node_color_ = color; }
void DagDrawingTool::label_color( int color )
{ label_color_ = color; }
void DagDrawingTool::background_color( int color )
{ background_color_ = color; }
int DagDrawingTool::primary_link_color() const
{ return prim_link_color_; }
int DagDrawingTool::secondary_link_color() const
{ return sec_link_color_; }
int DagDrawingTool::host_para_link_color() const
{ return host_para_color_; }
int DagDrawingTool::bad_link_color() const
{ return bad_link_color_; }
int DagDrawingTool::node_color() const
{ return node_color_; }
int DagDrawingTool::label_color() const
{ return label_color_; }
int DagDrawingTool::background_color() const
{ return background_color_; }
	
int DagDrawingTool::window_id() const
{ return window_id_; }
*/

/*
const char* DagDrawingTool::name( int row, int index )
{
	DagNodeTable& table = *node_table;
	return name( table[row][index] );
}
*/
/*
int DagDrawingTool::name_width( int row )
{
	DagNodeTable& table = *(node_table);
	ModelEntity* node_ptr = 0;
	if( table[row].length() > 0 )
		node_ptr = table[row][0];
	return name_width( node_ptr );
}
*/
	
/*
int DagDrawingTool::name_width( ModelEntity* node_ptr )
{
  if( node_ptr == NULL )
     return 4;
  else if( CAST_TO( node_ptr, Body ) ||
           CAST_TO( node_ptr, RefVolume ) ||
           CAST_TO( node_ptr, RefFace ) ||
           CAST_TO( node_ptr, RefEdge ) ||
           CAST_TO( node_ptr, RefVertex ) )
  {
    switch( ref_entity_label_type_ )
    {
      case DDT_LABEL_REFENTITY_WITH_NAME: 
         return DDT_DEFAULT_NAME_LEN;
      case DDT_LABEL_REFENTITY_WITH_TYPE: 
         return DDT_DEFAULT_TYPE_LEN + DDT_DEFAULT_ID_LEN;
      case DDT_LABEL_REFENTITY_WITH_ABBR: 
         return DDT_DEFAULT_ABBREV_LEN + DDT_DEFAULT_ID_LEN - 1;
      case DDT_LABEL_REFENTITY_WITH_ID: 
         return DDT_DEFAULT_ID_LEN;
      default: 
         assert( ! ref_entity_label_type_ );
    }
  }
	else
  {
    switch( other_entity_label_type_ )
    {
      case DDT_LABEL_OTHER_WITH_TYPE: return DDT_DEFAULT_TYPE_LEN;
      case DDT_LABEL_OTHER_WITH_ABBR: return DDT_DEFAULT_ABBREV_LEN;
      case DDT_LABEL_OTHER_NONE: return 1;
      default: assert( ! other_entity_label_type_ );
    }
	}
	
	return 0;
}
*/
/*
const char* DagDrawingTool::name( ModelEntity* node_ptr )
{
	static char buffer[32];
	
	if( ! node_ptr )
	{
			strcpy( buffer, "NULL" );
			return buffer;
	}
	
	ModelEntity* entity_ptr = node_ptr;
	if( ! entity_ptr )
	{
		strcpy( buffer, "NODE" );
		return buffer;
	}
	
	CubitString s;
	RefEntity* re_ptr = CAST_TO( entity_ptr, RefEntity );
	if( re_ptr )
	{
		switch( ref_entity_label_type_ )
		{
			case DDT_LABEL_REFENTITY_WITH_NAME:
				s = re_ptr->entity_name();
				strncpy( buffer, s.c_str(), 31 );
				break;
			case DDT_LABEL_REFENTITY_WITH_TYPE:
				sprintf( buffer, "%0.24s %d",
                                                 re_ptr->class_name(), re_ptr->id() );
				break;
			case DDT_LABEL_REFENTITY_WITH_ABBR:
				sprintf( buffer, "%0.24s%d", short_name( re_ptr ), re_ptr->id() );
				break;
			case DDT_LABEL_REFENTITY_WITH_ID:
				sprintf( buffer, "%d", re_ptr->id() );
				break;
			default: assert( ! ref_entity_label_type_ );
		}
	}
	else
	{
		switch( other_entity_label_type_ )
		{
			case DDT_LABEL_OTHER_WITH_TYPE:
				strncpy( buffer,  
                 entity_ptr->class_name( ), 31 );
				break;
			case DDT_LABEL_OTHER_WITH_ABBR:
				strncpy( buffer, short_name( entity_ptr ), 31 );
				break;
			case DDT_LABEL_OTHER_NONE:
				buffer[0] = '\0';
				break;
			default: assert( ! other_entity_label_type_ );
		}
	}
	return buffer;
}
*/



void DagDrawingTool::get_relatives( ModelEntity* source_ptr,
	DLIList<ModelEntity*>& result_set, int direction )
{
	result_set.clean_out();
	if( direction == DDT_UP_DAG ) 
		source_ptr->get_parents( &result_set );
	else if( direction == DDT_DOWN_DAG )
	  source_ptr->get_children( &result_set  );
	else assert( (direction == DDT_UP_DAG) || (direction == DDT_DOWN_DAG) );
}

void DagDrawingTool::get_relatives( DLIList<ModelEntity*>& source_set,
	DLIList<ModelEntity*>& result_set, int direction )
{
	DLIList<ModelEntity*> temp_set;
	result_set.clean_out();
	for( int i = 0; i < source_set.size();  i++ )
	{
		get_relatives( source_set.get_and_step(), temp_set, direction );
		result_set.merge_unique( temp_set );
	}
}


int DagDrawingTool::link_type( int source_row, int source_index, 
                               int target_row, int target_index )
{
	DagNodeTable& table = *node_table;
	if( abs( source_row - target_row ) != 1 ) return DDT_LINK_NONE;
	return link_type( table[source_row][source_index], 
	                  table[target_row][target_index], 
										source_row < target_row ? DDT_UP_DAG : DDT_DOWN_DAG );
}
int DagDrawingTool::link_type( ModelEntity* source_node,
	ModelEntity* target_node, int direction )
{
	int result = DDT_LINK_NONE;
  DLIList<ModelEntity*> parents, children;
	if( direction == DDT_UP_DAG )
	{
    source_node->get_parents(&parents);
    target_node->get_children(&children);
		bool down = parents.is_in_list(target_node);
    bool up   = children.is_in_list(source_node);
    if ( !down && !up )
      result = DDT_LINK_NONE;
    else if( down && up )
      result = DDT_LINK_BOTH;
    else
      result = DDT_BAD_LINK;
  }
	else if( direction == DDT_DOWN_DAG )
	{
    source_node->get_children(&children);
    target_node->get_parents(&parents);
		bool up   = children.is_in_list(target_node);
    bool down = parents.is_in_list(source_node);
    if ( !down && !up )
      result = DDT_LINK_NONE;
    else if( down && up )
      result = DDT_LINK_BOTH;
    else
      result = DDT_BAD_LINK;
	}
	return result;
}



DagNodeRow::DagNodeRow( DLIList<ModelEntity*>& node_list )
{
	length_ = node_list.size();
	array_ = 0;
	if( length_ > 0 )
	{
		array_ = new ModelEntity*[length_];
		node_list.reset();
		for( int i = 0; i < length_; i++ )
			array_[i] = node_list.get_and_step();
	}
}

DagNodeRow::~DagNodeRow()
{
	if( array_ ) delete [] array_;
	array_ = 0;
	length_ = 0;
}

ModelEntity*& DagNodeRow::operator[]( int index )
{
	assert( (index >= 0) && (index < length_) );
	return array_[index];
}





DagNodeTable::~DagNodeTable( )
{
	if( array_ )
	{
		for( int i = 0; i < length_; i++ )
		{
			if( array_[i] ) 
			{
				delete array_[i];
				array_[i] = 0;
			}
		}
		delete [] array_;
	}
	array_ = 0;
	length_ = 0;
}

DagNodeRow& DagNodeTable::operator[]( int row )
{
	assert( (row >= 0) && (row < length_) && array_[row] );
	return *(array_[row]);
}





class DLT_IdTable
{
public:
  DLT_IdTable( ){ current_insert_pos = 0; }
  ~DLT_IdTable( ){ }
  
  int find_id( ModelEntity* ME_ptr );
  
private:
  ModelEntity* ME_ptrs_[ 500 ];
  int current_insert_pos;
};

void DagDrawingTool::printDag(DLIList<RefFace*> &face_list, int depth) 
{
  DLIList<ModelEntity*> entity_list;
  CAST_LIST_TO_PARENT(face_list, entity_list);
  printDag(entity_list, depth);
}

void DagDrawingTool::printDag(DLIList<Body*> &body_list, int depth)
{
  DLIList<ModelEntity*> entity_list;
  CAST_LIST_TO_PARENT(body_list, entity_list);
  printDag(entity_list, depth);
}

void DagDrawingTool::printDag(DLIList<ModelEntity*> &entity_list, int depth) 
{
  int i;
  for (i = entity_list.size(); i > 0; i--) {
    printDag(entity_list.get_and_step(), depth);
  }
}

void DagDrawingTool::printDag( ModelEntity* ME_ptr, int depth )
{
  assert( ME_ptr != NULL );
  ModelEntity* node = ME_ptr;
  print_dag( node, depth );
}

void DagDrawingTool::print_dag( ModelEntity* any_dag_node, int depth )
{
	if( depth >= 0 ) print_node( any_dag_node, 0, "**" );
	print_dag( any_dag_node, 1, depth );
	if( depth < 0 ) print_node( any_dag_node, 0, "**" );
}

void DagDrawingTool::print_node( ModelEntity* node_ptr, int indent, const char* prefix ) 
{
	const char* name_ptr = typeid(*node_ptr).name();
  while (isdigit(*name_ptr)) name_ptr++;

  int id = get_id( node_ptr );
	PRINT_INFO("%*c%s %s %d", indent*3, ' ', prefix, name_ptr, id);
  
  PRINT_INFO("\n");
}

int DagDrawingTool::get_id( ModelEntity* node_ptr )
{
  static DLT_IdTable SE_table;
  static DLT_IdTable GpE_table;

	RefEntity*       re_ptr = CAST_TO( node_ptr, RefEntity );
	GroupingEntity* gpe_ptr = CAST_TO( node_ptr, GroupingEntity );
	SenseEntity*     se_ptr = CAST_TO( node_ptr, SenseEntity );
	if( re_ptr )      return re_ptr->id();
	else if( gpe_ptr ) return GpE_table.find_id( gpe_ptr );
	else if( se_ptr )  return SE_table.find_id( se_ptr );
	else               return 0;
}

void DagDrawingTool::print_dag( ModelEntity* node_ptr, int indent, int depth )
{
	DLIList<ModelEntity*> relatives;
	if( depth > 0 ) {
		node_ptr->get_children( &relatives );
		for( int i = relatives.size(); i > 0; i-- ) {
			ModelEntity* child_ptr = relatives.get_and_step();
				
			print_node( child_ptr, indent, "<>" );
			print_dag( child_ptr, indent + 1, depth - 1 );
		}
	}
	else if( depth < 0 ) {
		node_ptr->get_parents( &relatives );
		for( int i = relatives.size(); i > 0; i-- ) {
			ModelEntity* parent_ptr = relatives.get_and_step();
				
			print_dag( parent_ptr, indent + 1, depth + 1 );
			print_node( parent_ptr, indent, "<>" );
		}
	}
}

int DLT_IdTable::find_id( ModelEntity* ME_ptr )
{
  ME_ptr = CAST_TO( ME_ptr, ModelEntity );
  if( ME_ptr == NULL ) return -1;
  
  for( int i = 0; i < current_insert_pos; i++ )
  {
    if( ME_ptrs_[i] == ME_ptr ) return i + 1;
  }
  if( current_insert_pos == 500 ) return -1;
  ME_ptrs_[ current_insert_pos ] = ME_ptr;
  current_insert_pos++;
  return current_insert_pos;
}



