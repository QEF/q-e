#!/usr/bin/awk -f

# prediff.awk -- preprocessor for checking outputs of PWscf examples
# checking is done in three steps: preprocess, diff against reference
# data, postprocess
# this way "false positives" are eliminated, that is, differences that
# don't mean that something went wrong

# for each (group of) line(s)
{ check_line(); }

function check_line()
{
  # mark (groups of) lines that may be "false positives", by prepending
  # a key of the form "@key@"
  # postprocessor will check them based on key

  if (match($0, "Today is") || match($0, "cpu time") || match($0, "CPU"))
    {
      print_key("TIMING");
      if (getline)
	{
	  if (NF == 0)
	    print_key("TIMING");
	  else
	    check_line();
	}
    }
  else if (match($0, "Parallel version"))
    {
      print_key_n("PARALLEL", 2);

      # see /Reading header from file/ below
      getline; $0 = ""; print_key("PARALLEL");

      while (getline && NF == 0)
	print_key("PARALLEL");
      print;
    }
  else if (match($0, "Planes per process"))
    {
      print_key_n("PARALLEL", 2);
    }
  else if (match($0, "Proc/  planes cols"))
    {
      print_key("PARALLEL");
      while (getline && NF > 0)
	print_key("PARALLEL");
      print_key("PARALLEL");

      # see /Reading header from file/ below
      getline; $0 = ""; print_key("PARALLEL");
    }
  else if (match($0, "Parallel routines"))
    {
      print_key("PARALLEL");
    }
  else if (match($0, "npwx   =") || match($0, "ngl    ="))
    {
      print_key("PARALLEL");
    }
  else if (match($0, "Reading header from file") \
	   || match($0, "Reading data from file"))
    {
      # must erase these lines because of a subtle problem with
      # blank lines before/after parallel data
      $0 = ""; print;
    }
  else if (match($0, "FFT grid") || match($0, "smooth grid"))
    {
      print_key("FFTGRID");
    }
  else if (match($0, "wk ="))
    {
      print_key("KWEIGHT");
    }
  else if (match($0, "Self-consistent Calculation"))
    {
      print_key("ITERATION");
      while (getline && ! match($0, "End of self-consistent calculation"))
	print_key("ITERATION");
      print_key("ITERATION");
      print "@CHECKPOINT@";
    }
  else if (match($0, "Band Structure Calculation"))
    {
      print_key("ITERATION");
      while (getline && ! match($0, "End of band structure calculation"))
	print_key("ITERATION");
      print_key("ITERATION");
    }
  else if (match($0, "band energies") || match($0, "bands"))
    {
      print;
      getline; print;
      while (getline && NF > 0)
	print_key("BANDS");
      print;
    }
  else if (match($0, "the Fermi energy is"))
    {
      print_key("FERMI");
    }
  else if (match($0, "total energy"))
    {
      print_key_n("ENERGY", 2);
    }
  else if (match($0, "band energy sum") || match($0, "contribution") \
	   || match($0, "correction for metals"))
    {
      print_key("CONTRIBUTION");
    }
  else if (match($0, "total   stress"))
    {
      print_key("PRESSURE");
      getline; print_key_n("STRESS", 3);
    }
  else if (match($0, "ddv_scf"))
    {
      print_key("DDV");
    }
  else if (match($0, "Dielectric constant"))
    {
      print;
      getline; print;
      getline; print_key_n("DIELECTRIC", 3);
    }
  else if (match($0, "omega"))
    {
      print_key("OMEGA");
    }
  else if (match($0, "atom.*type.*force ="))
    {
      print_key("FORCE");
    }
  else if (match($0, "Total force ="))
    {
      print_key("TFORCE");
    }
  else if (match($0, "energy new *="))
    {
      print_key("ITERATION");
      while (getline && ! match($0, "bfgs converged"))
	print_key("ITERATION");
      print_key("ITERATION");
      print "@CHECKPOINT@";
    }
  else if (match($0, "Final energy *="))
    {
      print_key("ENERGY");
    }
  else if (match($0, "Ekin =.*T =.*Etot ="))
    {
      print_key("TEMPERATURE");
    }
  else if (match($0, "Search of equilibrium positions:"))
    {
      print_key("ITERATION");
      while (getline && ! match($0, "convergence achieved, Efinal="))
	print_key("ITERATION");
      print "@CHECKPOINT@";
      print_key("EFINAL");
    }
  else if (match($0, "new positions"))
    {
      while (NF > 0)
	{
	  print;
	  while (getline && NF == 4)
	    print_key("POSITIONS");
	}
      print;

      # geometry relaxation starting
      while (getline && ! match($0, "convergence achieved, Efinal="))
	print_key("ITERATION");
      print "@CHECKPOINT@";
      print_key("EFINAL");
    }
  else if (match($0, "ATOMIC_POSITIONS"))
    {
      print;
      while (getline && NF == 4)
	print_key("POSITIONS");
      print;
    }
  else
    {
      # unrecognized type of line, print as is
      print;
    }
}

function print_key(key)
{
  print "@" key "@ " $0;
}

function print_key_n(key, n)
{
  print_key(key);
  for (i=1; i<n; i++)
    {
      getline; print_key(key);
    }
}
