/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ReweightBase.h"
#include "tools/File.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_ITS
/*
Calculate the weights configurations should contribute to the histogram in a simulation in which a ITS enhanced sampling acts upon the system.

This command allows you to use the reweighting algorithm discussed in \cite PratyushReweighting when constructing a histogram
of the configurations visited during a metadynamics simulation.

\par Examples

In the following example there is a metadynamics bias acting on the distance between atoms 1 and 2.  Clearly, this 
bias will have an effect on the region of phase space that will be sampled when an MD simulation is
run using this variable.  Consequently, when the histogram as a function of the angle, \f$a\f$, is accumulated,
we use reweighting into order to discount the effect of the bias from our final histogram.  We do not use 
\ref REWEIGHT_BIAS here, however, as the bias changes with time.  We thus use the reweighting algorithm for 
metadynamics instead.  Notice also that we have to specify how often we would like to calculate the c(t) reweighting
factor and the grid over which we calculate c(t) in the input to the METAD command.

\verbatim
a: ANGLE ATOMS=1,2,3
x: DISTANCE ATOMS=1,2
METAD ARG=x PACE=100 SIGMA=0.1 HEIGHT=1.5 BIASFACTOR=5 GRID_MIN=0 GRID_MAX=10 GRID_BIN=100 REWEIGHTING_NGRID=100 REWEIGHTING_NHILLS=50 

as: REWEIGHT_ITS SIM_TEMP=300

HISTOGRAM ...
  ARG=a
  GRID_MIN=0.0
  GRID_MAX=pi
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
  LABEL=hB
... HISTOGRAM

DUMPGRID GRID=hB FILE=histoB STRIDE=1 FMT=%8.4f

\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightITS : public ReweightBase {
private:
    unsigned narg;
    unsigned effid;
    unsigned potid;
/// The temperature at which you are running the simulation
	bool use_fb;
	bool calc_rct;

	double kB;
	double sim_temp;
	double rw_temp;
	double beta0;
	double betaRW;
	double shift_const;
	double rct;
	
	std::vector<double> int_temps;
	std::vector<double> betak;
	std::vector<double> fb;
	
	std::string fb_file;
	
	//~ unsigned nreplica;
	
	unsigned read_fb_file(const std::string& fname);
	double calc_rct_value(double rwtemp);
public:
  static void registerKeywords(Keywords&);
  explicit ReweightITS(const ActionOptions&ao);
  double getLogWeight() const ;
};

PLUMED_REGISTER_ACTION(ReweightITS,"REWEIGHT_ITS")

void ReweightITS::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","The arguments"); 
  keys.add("optional","RW_TEMP","the reweight temperature. If it is not assigned, it will be equal to the system temperatures");
  keys.add("optional","SHIFT","a constant energy to shift the bias potential when they are too large or too small");
  keys.add("optional","RCT","the c(t) at the reweight temperature");
}

ReweightITS::ReweightITS(const ActionOptions&ao):
Action(ao),
ReweightBase(ao),
narg(getNumberOfArguments()),
shift_const(0),
rct(0)
{
	plumed_massert(narg>=2,"The reweight of ITS need at least two arguments: the potential energy of the system and the effective potential\n");
	kB=plumed.getAtoms().getKBoltzmann();
	
	beta0=1.0/simtemp;
	sim_temp=simtemp/kB;
	
	rw_temp=-1;
	parse("RW_TEMP",rw_temp);
	if(rw_temp>0)
		betaRW=1.0/(kB*rw_temp);
	else
	{
		betaRW=beta0;
		rw_temp=sim_temp;
	}

	parse("RCT",rct);
	parse("SHIFT",shift_const);
	
	log.printf("  with Boltzmann constant: %f\n",kB);
	log.printf("  with simulation temperature: %f K (%f)\n",sim_temp,beta0);
	log.printf("  with reweight temperature: %f K (%f)\n",rw_temp,betaRW);
	log.printf("  with revise factor c(t): %f\n",rct);
	log.printf("  with constant shift energy: %f\n",shift_const);
}

double ReweightITS::getLogWeight() const
{
	double pot=getArgument(0);
	double eff=0;
	for(unsigned i=1;i<narg;++i)
		eff+=getArgument(i);

	return beta0 * eff + (shift_const - pot - rct ) * betaRW;
}

}
}
