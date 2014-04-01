//- Class: IGUIObservers
//- Description: interface for GUI notification from CubitEntities
//- Created: 06/18/02 by Byron Hanks, Sandia National Laboratories
//- Checked By:
//- Version:
class IGUIObservers
{
public:
    // Notify the GUI that the entity holding the IGUIObservers pointer
    // has been destroyed.
    virtual void Destroyed() = 0;
};

