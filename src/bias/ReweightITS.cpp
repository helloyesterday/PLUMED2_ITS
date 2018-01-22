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
/// The temperature at which you are running the simulation
	bool rw_other;
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
	
	unsigned nreplica;
	
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
  keys.remove("TEMP");
  keys.add("compulsory","ARG","the biases that must be taken into account when reweighting. If you want to reweight the system at other temperatures, please set the ITS-energy at the 1st argument, the ITS-bias at the 2nd argument"); 
  keys.add("compulsory","SIM_TEMP","the reweighting temperature.  If it is not assigned, it will be equal to the system temperatures");
  keys.add("optional","RW_TEMP","the simulation temperature.");
  keys.add("optional","FB_FILE","using fb file to calculate reweighting factor"); 
  keys.addFlag("USE_FB_FILE",false,"use fb file to calculate the reweight. The eneryg of ITS must be the first term at ARG and DO NOT put the bias of ITS at ARG");
  keys.addFlag("CALC_RCT",false,"calculate the c(t) from the fb file and use it at the reweighting");
  keys.add("optional","SHIFT","a constant energy to shift the bias potential when they are too large or too small");
}

ReweightITS::ReweightITS(const ActionOptions&ao):
Action(ao),
ReweightBase(ao),
rw_other(false),
shift_const(0),
rct(0)
{
	
	kB=plumed.getAtoms().getKBoltzmann();
	parse("SIM_TEMP",sim_temp);
	
	beta0=1.0/(kB*sim_temp);
	
	rw_temp=-1;
	parse("RW_TEMP",rw_temp);
	if(rw_temp>0)
	{
		betaRW=1.0/(kB*rw_temp);
		rw_other=true;
	}
	else
	{
		betaRW=beta0;
		rw_temp=sim_temp;
	}
	
	parse("SHIFT_CONST",shift_const);
	
	plumed_massert(getNumberOfArguments()>=2,"to reweight at the other temperature need at least two arguments: its.energy and its.bias. Please put them as the 1st and 2nd ARG");
	
	parse("FB_FILE",fb_file);
	parseFlag("USE_FB_FILE",use_fb);
	if(use_fb&&fb_file.size()==0)
		plumed_merror("USE_FB_FILE must set the fb file at FB_FILE");
	parseFlag("CALC_RCT",calc_rct);
	if(calc_rct&&fb_file.size()==0)
		plumed_merror("CALC_RCT must set the fb file at FB_FILE");
	
	if(use_fb||calc_rct)
		read_fb_file(fb_file);
	if(calc_rct)
		rct=calc_rct_value(rw_temp);
}

double ReweightITS::getLogWeight() const
{
	if(use_fb)
	{
		double its_energy=getArgument(0);
		double other_bias=shift_const;
		for(unsigned i=1;i<getNumberOfArguments();++i)
			other_bias+=getArgument(i);
	
		double factor=0;
		for(unsigned i=0;i!=nreplica;++i)
			factor+=exp(fb[i]+(betaRW-betak[i])*its_energy);
		return -std::log(factor) + betaRW * other_bias - rct;
	}
	else
	{
		if(rw_other)
		{
			double its_energy=getArgument(0);
			double its_bias=getArgument(1);
			double other_bias=shift_const;
			for(unsigned i=2;i<getNumberOfArguments();++i)
				other_bias+=getArgument(i);

			return beta0 * its_bias + ( beta0 - betaRW ) * its_energy + betaRW * other_bias - rct;
		}
		else
		{
			double bias=shift_const;
			for(unsigned i=0;i<getNumberOfArguments();++i)
				bias+=getArgument(i);
			return betaRW * bias - rct;
		}
	}
}

unsigned ReweightITS::read_fb_file(const std::string& fname)
{	
	IFile ifb;
		
	if(!ifb.FileExist(fname))
		plumed_merror("Cannot find fb file " + fname );
    ifb.open(fname);
    
	ifb.allowIgnoredFields();

	unsigned read_count=0;
	bool is_set_ratios=false;
	
	std::vector<double> int_ratios,int_temps;
	while(ifb)
	{
		
		std::string coe_type="TEMPERATURE";
		if(ifb.FieldExist("COEFFICIENT_TYPE"))
			ifb.scanField("COEFFICIENT_TYPE",coe_type);

		if(coe_type=="RATIO")
			is_set_ratios=true;
		else if(coe_type!="TEMPERATURE")
			plumed_merror("unrecognized COEFFICIENT_TYPE "+coe_type);

		
		if(ifb.FieldExist("NREPLICA"))
		{
			int int_ntemp;
			ifb.scanField("NREPLICA",int_ntemp);
			nreplica=int_ntemp;
		}
		else
			plumed_merror("Can't find the NREPLICA from the fb file: "+fname);
		
		if(is_set_ratios)
			int_ratios.resize(nreplica);
		else
			int_temps.resize(nreplica);
		fb.resize(nreplica);

		if(ifb)
		{
			if(!ifb.FieldExist("fb_value"))
				plumed_merror("cannot found \"fb_value\" in file \""+fname+"\"");
			if(!is_set_ratios&&!ifb.FieldExist("temperature"))
				plumed_merror("cannot found \"temperature\" in file \""+fname+"\"");
			if(is_set_ratios&&!ifb.FieldExist("replica_ratio"))
				plumed_merror("cannot found \"replica_ratio\" in file \""+fname+"\"");
		}
		else
			break;
		
		double tmpfb;
		for(unsigned i=0;i!=nreplica;++i)
		{
			ifb.scanField("fb_value",tmpfb);
			double tmpt;
			if(is_set_ratios)
			{
				ifb.scanField("replica_ratio",tmpt);
				int_ratios[i]=tmpt;
				int_temps[i]=sim_temp/tmpt;
			}
			else
			{
				ifb.scanField("temperature",tmpt);
				int_temps[i]=tmpt;
				int_ratios[i]=sim_temp/tmpt;
			}
			
			fb[i]=tmpfb;
			ifb.scanField();
		}
		++read_count;
	}
	
	betak.resize(nreplica);
	if(is_set_ratios)
		int_temps.resize(nreplica);
	for(unsigned i=0;i!=nreplica;++i)
	{
		if(is_set_ratios)
		{
			betak[i]=beta0*int_ratios[i];
			int_temps[i]=sim_temp/int_ratios[i];
		}
		else
			betak[i]=1.0/(kB*int_temps[i]);
	}

	ifb.close();
	return read_count;
}

double ReweightITS::calc_rct_value(double rwtemp)
{
	double dtl=0,dth=0;
	
	unsigned rwid=0;
	for(unsigned i=0;i!=nreplica-1;++i)
	{
		dtl=rwtemp-int_temps[i];
		dth=int_temps[i+1]-rwtemp;
		if(dtl>=0&&dth>=0)
		{
			rwid=i;
			break;
		}
	}
	
	if(rwid+1==nreplica)
		plumed_merror("Can't find the fb value responds to the reweight temperature");
		
	double fbl=fb[rwid];
	double fbh=fb[rwid+1];
	double dfb=fbh-fbl;
	
	double rwfb;
	if(dtl<dth)
		rwfb=fbl+dfb*dtl/(dtl+dth);
	else
		rwfb=fbh-dfb*dth/(dtl+dth);
		
	return -(rwfb+std::log(double(nreplica)))/betaRW;;
}

}
}
