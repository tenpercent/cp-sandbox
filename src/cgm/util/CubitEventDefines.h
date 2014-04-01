#ifndef __CUBIT_EVENT_DEFINES__
#define __CUBIT_EVENT_DEFINES__

/* Event types. These are used in the Subject-Observer patterns in the
 * code to specify what event occurred.  This allows the Observers to
 *          update themselves accordingly, if required. */
enum CubitEventType { INVALID_EVENT_TYPE = -1,
                 VGI_BODY_CONSTRUCTED,
                   /* A ModelEntity was created */
                 MODEL_ENTITY_CONSTRUCTED,
                   /* A ModelEntity was created */
                 MODEL_ENTITY_DESTRUCTED,
                   /* A ModelEntity was deleted */
                 COMPARISON_FOUND,
                   /* A partner entity was found (they compare successfully) */
                 ENTITY_SURVIVED_MERGE,
                   /* An entity was involved in a merge and survived */
                 MERGE_COMPLETED,
                   /* A merge operation is completed successfully */
                 MERGE_ABORTED,
                   /* A merge operation was aborted */
                 MESH_DELETED,
                   /* Mesh associated with 1 or more RefEntities was deleted */
                 MESH_CREATED,
                   /* 1 or more RefEntities was meshed */
                 MESH_MODIFIED,
                   /* 1 or more RefEntities had its exiting mesh modified */
                 GEOMETRY_TOPOLOGY_MODIFIED,
                   /* Both geometry and topology was
                    modified. e.g. partitioned.
                   */
                 TOPOLOGY_MODIFIED,
                   /* The topology of a RefVolume or one of its sub-entities
                    * was changed.  Used for virtual topology. */
                 GEOMETRY_MODIFIED,
                   /* The geometry of a RefEntity was altered. */
                 HEALER_COMPLETED,
                   /* Notifies that the healer completed */
                 NEW_ENTITY_UNMERGED,
                   /*A surface, curve, or vertex was unmerged. */
                 DAG_NODE_DESTRUCTED,
                   /* A DAG Node was destructed
                    * CHILD_TO_BE_SWAPPED
                    * A child is about to be swapped */
                 FREE_REF_ENTITY_GENERATED,
                   /* A vertex that is not a part of a curve,
                    * or a curve that is not a part of a surface,
                    * was just generated */
                 ENTITY_CONSTRUCTED,
                 ENTITY_DESTRUCTED,
                   /* An generic entity is constructed or destructed; used to
                    * notify observers of destruction of observables */
                 TOP_LEVEL_ENTITY_DESTRUCTED,
                   /* A body or free entity and all it's children
                    * are about to be destructed.
                    * The pointer is still valid at this point for all
                    * calls except geometry queries. 
                    * The pointer can be cast to a higher level entity
                    * if desired. */

                 ENTITY_NAME_CHANGED,
                   /* the name for an entity changed */
                 
                 ENTITY_GEOMETRY_COLOR_CHANGED,
                 ENTITY_MESH_COLOR_CHANGED,
                   /* the color for an entity changed */
                 ENTITY_VISIBILITY_CHANGED,
                  /* the visibility of an entity changed */

                 ENTITIES_MERGED,
                   /* two entities are merged together. See MergeEvent. */
                 ID_SET,
                   /* an entity had its id set (changed). See IdSetEvent. */
                 IDS_COMPRESSED,
                   /* ids have been compressed. See model.cpp */
                 MODEL_ENTITY_HIDDEN,
                   /* A sub-type of MODEL_ENTITY_DESTRUCTED to indicate
                    * that the entity in question is being hidden (VG)
                    * rather than actually being destructed.           */
                 MODEL_ENTITY_RESTORED,
                   /* The reverse of MODEL_ENTITY_HIDDEN, and a sub-type
                      of MODEL_ENTITY_CONSTRUCTED */
                 MODEL_RESET,
                   /* Indicates a 'reset' command has been issued and commpleted. */
                 GROUP_MODIFIED,
                   // Group was modified.  Sent only once after all
                   // modifications are complete.
                 GRAPHICS_MODE_MODIFIED,
                   // The overall graphics mode has been changed in some way.
                   // This could result from the perspective being changed, the
                   // rendering mode changed, the pick mode changed, etc.
                 
                 GENESIS_ENTITY_CREATED,
                   // Block, nodeset or sideset was created.
                 GENESIS_ENTITY_DELETED,
                   // Block, nodeset or sideset was deleted.
                 GENESIS_ENTITY_MODIFIED,
                   // Block, nodeset or sideset was modified.
                 SUSPEND_GENESIS_PROCESSING,
                   // Suspend processing (in CubitInterface) for Genesis events
                 RESUME_GENESIS_PROCESSING,
                   // Resume processing (in CubitInterface) for Genesis events
                 UPDATE_GENESIS_DISPLAY,
                   // Force an update of all genesis entity display


                   // ******** Assembly Events  *********
                   //   All assembly events pass an AssemblyEvent
                   //   object.  The functions used to get valid
                   //   data for a particular event type is listed
                   //   right after the event type name.
                 
                 ASSEMBLY_ADD_CHILD,
                   // (get_assembly(), get_node())
                   // Sent just after a new child is added to an assembly
                 
                 ASSEMBLY_REMOVE_CHILD,
                   // (get_assembly(), get_node())
                   // Sent just after a child is removed from an assembly
                 
                 ASSEMBLY_NODE_NAME_CHANGE,
                   // (get_node())
                   // Sent just after an assembly's name is changed. Note
                   // that the instance number may change along with the
                   // name, but this won't trigger a separate
                   // ASSEMBLY_NODE_INSTANCE_CHANGE event.
                 
                 ASSEMBLY_NODE_INSTANCE_CHANGE,
                   // (get_node())
                   // Sent just after an assembly or part's instance
                   // number is changed
                 
                 ASSEMBLY_NODE_STRING_PROPERTY_CHANGE,
                   // (get_node(), get_property())
                   // Sent just after an assembly or node has one of its
                   // properties changed, and that property is a string
                   // value
                 
                 PART_ADD_VOLUME,
                   // (get_part(), get_volume())
                   // Sent just after a volume is associated with a part
                 
                 PART_REMOVE_VOLUME,
                   // (get_part(), get_volume())
                   // Sent just after a volume is removed from a part
                 
                 ASSEMBLY_TREE_STRING_PROPERTY_CHANGE
                   // (get_property())
                   // Sent just after a property of the tree changes.
                   // These are properties that apply to the tree as a whole,
                   // not to any individual AssemblyNode.
                 
                   // 
                   // ********* End of Assembly Events **********
};

#endif // __CUBIT_EVENT_DEFINES__
