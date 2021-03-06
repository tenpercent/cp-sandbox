
      double precision eta_cutoff, etamax
      double precision thmin, thmax

      common /nitparam/ eta_cutoff, etamax, thmin, thmax
      save /nitparam/

c
c nitparam contains some parameters that control the nonlinear
c iterations.  In some cases, the default values reflect prevailing  
c practice; in other cases, they are chosen to produce good 
c average-case behavior.  To change the default values, include this 
c common block in the main program and set the desired variables 
c according to the following:
c
c    eta_cutoff   - parameter used to determine when to disable 
c                   safeguarding the update of the forcing term.  It
c                   only has meaning when ieta .ne. 3.  The default
c                   value is 0.1.  A value of 0.0 will enable 
c		    safeguarding always; a value of 1.0 will disable 
c		    safeguarding always. 
c
c    etamax       - parameter used to provide an upper bound on the 
c		    forcing terms when input(10) .ne. 3. This is 
c		    necessary to ensure convergence of the inexact Newton 
c		    iterates and is imposed whenever eta would otherwise 
c		    be too large. (An overly large eta can result from 
c		    the updating formulas when input(10) .ne. 3 or from 
c                   safeguarding when the previous forcing term has been 
c		    excessively increased during backtracking.) The 
c		    default value of etamax is 1.0 - 1.e-4.  When 
c		    backtracking occurs several times during a nonlinear 
c		    solve the forcing term can remain near etamax for several
c                   nonlinear steps and cause the nonlinear iterations
c                   to nearly stagnate.  In such cases a smaller value of 
c                   etamax may prevent this.  Valid values are in the 
c                   range (0.0, 1.0).
c
c    thmin        - when backtracking occurs, this is the smallest
c                   reduction factor that will be applied to the current
c                   step in a single backtracking reduction.  The default
c                   value is 0.1.  Valid  values are in the range
c                   [0.0, thmax].
c
c    thmax        - when backtracking occurs, this is the largest
c                   reduction factor that will be applied to the current
c                   step in a single backtracking reduction.  The default
c                   value is 0.5.  Valid values are in the range
c                   [thmin, 1.0).
c
