/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#include "GetLongOpt.hpp"
#include "CubitMessage.hpp"
#include "CubitUtil.hpp"
const char* const GetLongOpt::TOGGLE_ON = "True";
const char* const GetLongOpt::TOGGLE_OFF = "False";

GetLongOpt::GetLongOpt(const char optmark)
{
   table = last = 0;
   ustring = "[valid options and arguments]";
   enroll_done = 0;
   optmarker = optmark;
}

GetLongOpt::~GetLongOpt()
{
   Cell *t = table;

   while ( t ) {
      Cell *tmp = t;
      t = t->next;
      if (tmp->value != NULL && 
          tmp->value != (char*) ~0 &&
          tmp->value != TOGGLE_ON &&
          tmp->value != TOGGLE_OFF)
        CubitUtil::util_strdup_free( const_cast<char*>(tmp->value) );
      delete tmp;
   }
}

const char *GetLongOpt::basename(const char *prog_name) const
{
   const char *s;

   s = strrchr(prog_name, '/');
   if ( s == NULL )
     s = prog_name;
   else ++s;

   return s;
}

int
GetLongOpt::enroll(const char * const opt, const OptType t,
const char * const desc, const char * const val)
{
   if ( enroll_done ) return 0;

   Cell *c = new Cell;
   c->option = opt;
   c->type = t;
   c->description = desc ? desc : "no description available";
   c->value = val ? CubitUtil::util_strdup(val) : 0;
   c->wasSet = 0;
   c->next = 0;

   if ( last == 0 ) {
      table = last = c;
   }
   else {
      last->next = c;
      last = c;
   }

   return 1;
}

const char *
GetLongOpt::retrieve(const char * const opt) const
{
   Cell *t;
   for ( t = table; t != 0; t = t->next ) {
      if ( strcmp(opt, t->option) == 0 )
	 return t->value;
   }
   PRINT_ERROR("GetLongOpt::retrieve - unenrolled option %c%s\n",
		   optmarker, opt );
   return 0;
}

int GetLongOpt::parse(int argc, char *const* argv)
{
  int this_optind = 1;
  
  pname = basename(*argv);
  enroll_done = 1;
  if ( argc-- <= 1 ) return this_optind;
  
  while ( argc >= 1 )
  {
    char *token = *++argv; --argc;
    
    if ( token[0] != optmarker || token[1] == optmarker )
      break;	/* end of options */
    
    ++this_optind;
    const char *tmptoken = ++token;
    while ( *tmptoken && *tmptoken != '=' )
      ++tmptoken;
      /* (tmptoken - token) is now equal to the command line option
         length. */
    
    size_t len = tmptoken - token;
    Cell *t;
    enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
    Cell *pc = 0;	// pointer to the partially-matched cell
    bool toggle = true; // toggle state for partial match cell
    for ( t = table; t != 0; t = t->next )
    {
      bool no = false;
      if ( (strncmp(t->option,token, len) == 0) ||
           (no = (t->type == Toggle && 
                  token[0]=='n'  && token[1] == 'o' && token[2] != '\0' && 
                  strncmp(t->option, token+2, len-2) == 0) ) ) 
      {
        if ( strlen(t->option) == (no ? len-2 : len) )
        {
            /* an exact match found */
          if (t->type == Toggle && *tmptoken == '\0' )
            tmptoken = no ? GetLongOpt::TOGGLE_OFF : GetLongOpt::TOGGLE_ON;
            
          int stat = setcell(t, tmptoken, *(argv+1), pname);
          if ( stat == -1 ) return -1;
          else if ( stat == 1 )
          {
            ++argv; --argc; ++this_optind;
          }
          matchStatus = ExactMatch;
          break;
        }
        else
        {
            /* partial match found */
          matchStatus = PartialMatch;
          pc = t;
          toggle = !no;
        }
      } /* end if */
    } /* end for */
    
    if ( matchStatus == PartialMatch )
    {
      if (pc->type == Toggle && *tmptoken == '\0' )
        tmptoken = toggle ? GetLongOpt::TOGGLE_ON : GetLongOpt::TOGGLE_OFF;
        
      int stat = setcell(pc, tmptoken, *(argv+1), pname);
      if ( stat == -1 ) return -1;
      else if ( stat == 1 ) 
      {
        ++argv; --argc; ++this_optind;
      }
    }
    else if ( matchStatus == NoMatch )
    {
      PRINT_ERROR("%s: unrecognized option %c%s\n",
                  pname, optmarker, strtok(token,"= \n") );
      return -1;		/* no match */
    }
    
  } /* end while */
  
  return this_optind;
}

int
GetLongOpt::parse(const char * const str_in, const char * p)
{
   enroll_done = 1;
   char* str = CubitUtil::util_strdup( str_in );
   char *token = strtok(str, " \t");
   const char *name = p ? p : "GetLongOpt";

   while ( token ) {
      if ( token[0] != optmarker || token[1] == optmarker ) {
        PRINT_ERROR("%s: nonoptions not allowed\n", name );
        CubitUtil::util_strdup_free(str);
        return -1;	/* end of options */
      }

      char *ladtoken = 0;	/* lookahead token */
      const char *tmptoken = ++token;
      while ( *tmptoken && *tmptoken != '=' )
        ++tmptoken;
      /* (tmptoken - token) is now equal to the command line option
        length. */

      size_t len = tmptoken - token;
      Cell *t;
      enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
      Cell *pc =0;	// pointer to the partially-matched cell
      bool toggle = true; // toggle state for partial match cell
      for ( t = table; t != 0; t = t->next ) {
        bool no = false;
      
        if ( (strncmp(t->option, token, len) == 0) ||
             (no = (t->type == Toggle && 
                    token[0]=='n'  && token[1] == 'o' && token[2] != '\0' && 
                    strncmp(t->option, token+2, len-2) == 0) ) ) 
        {
          if ( strlen(t->option) == (no ? len-2 : len) )
          {
              /* an exact match found */
            if (t->type == Toggle && *tmptoken == '\0' )
              tmptoken = no ? GetLongOpt::TOGGLE_OFF : GetLongOpt::TOGGLE_ON;

            ladtoken = strtok(0, " \t");
            int stat = setcell(t, tmptoken, ladtoken, name);
            if ( stat == -1 ) 
            {
              CubitUtil::util_strdup_free(str);
              return -1;
            }
            else if ( stat == 1 ) {
              ladtoken = 0;
            }
            matchStatus = ExactMatch;
            break;
          }
          else {
            /* partial match found */
            matchStatus = PartialMatch;
            pc = t;
            toggle = !no;
          }
        } /* end if */
      } /* end for */

      if ( matchStatus == PartialMatch ) {
        if (pc->type == Toggle && *tmptoken == '\0' )
          tmptoken = toggle ? GetLongOpt::TOGGLE_ON : GetLongOpt::TOGGLE_OFF;
        
        ladtoken = strtok(0, " \t");
        int stat = setcell(pc, tmptoken, ladtoken, name);
        if ( stat == -1 ) return -1;
        else if ( stat == 1 ) {
          ladtoken = 0;
        }
      }
      else if ( matchStatus == NoMatch ) {
        PRINT_ERROR("%s: unrecognized option %c%s\n",
        name, optmarker, strtok(token,"= \n"));
        CubitUtil::util_strdup_free(str);
        return -1;		/* no match */
      }

      token = ladtoken ? ladtoken : strtok(0, " \t");
   } /* end while */

   CubitUtil::util_strdup_free(str);
   return 1;
}

/* ----------------------------------------------------------------
GetLongOpt::setcell returns
   -1	if there was an error
    0	if the nexttoken was not consumed
    1	if the nexttoken was consumed
------------------------------------------------------------------- */
int
GetLongOpt::setcell(Cell *c, const char *valtoken, const char *nexttoken, 
                    const char *name)
{
   if ( c == 0 ) return -1;

   switch ( c->type ) {
    case GetLongOpt::Toggle :
      if ( valtoken != GetLongOpt::TOGGLE_ON && valtoken != GetLongOpt::TOGGLE_OFF ) {
        PRINT_ERROR("%s: unsolicited value for flag %c[no]%s\n",
			    name, optmarker, c->option );
        return -1;	/* unsolicited value specification */
      }
      c->value = valtoken;
      return 0;
      
    case GetLongOpt::Valueless :
      if ( *valtoken == '=' ) {
        PRINT_ERROR("%s: unsolicited value for flag %c%s\n",
			    name, optmarker,c->option );
        return -1;	/* unsolicited value specification */
      }
      if (!c->wasSet) {
        c->value = (c->value) ? 0 : (char *) ~0;
        c->wasSet = 1;
      }
      return 0;
    case GetLongOpt::OptionalValue :
      if ( *valtoken == '=' ) {
        c->value = CubitUtil::util_strdup(++valtoken);
      }
      else {
        if ( nexttoken != 0 && nexttoken[0] != optmarker ) {
          c->value = CubitUtil::util_strdup(nexttoken);
          return 1;
    	  }
      }

      // explicit return here, just to make sure another if-case isn't 
      // put in which falls through to the next case (in the absence of
      // a break statement)
      return 0;
    case GetLongOpt::MandatoryValue :
      int return_val;
      if ( *valtoken == '=' ) {
        c->value = CubitUtil::util_strdup(++valtoken);
        return_val = 0;
      }
      else {
        if ( nexttoken != 0 && nexttoken[0] != optmarker ) {
          c->value = CubitUtil::util_strdup(nexttoken);
          return_val = 1;
      }
      else {
        PRINT_ERROR("%s: mandatory value for %c%s\n",
          name, optmarker, c->option );
        return_val = -1;	/* mandatory value not specified */
        }
      }
      return return_val;
    default :
      break;
   }
   return -1;
}

void 
GetLongOpt::options(std::ostream &outfile) const
{
   Cell *t;

   for ( t = table; t != 0; t = t->next ) {
      outfile << "\t" << optmarker;
      if ( t->type == GetLongOpt::Toggle )
        outfile << "[no]";
      outfile << t->option;
      if ( t->type == GetLongOpt::MandatoryValue )
	 outfile << " <$val>";
      else if ( t->type == GetLongOpt::OptionalValue )
	 outfile << " [$val]";
      outfile << " (" << t->description << ")\n";
   }
}

void
GetLongOpt::usage( std::ostream &outfile) const
{
   outfile << "usage: " << pname << " " << ustring << "\n";
   options(outfile);
}

