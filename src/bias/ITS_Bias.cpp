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
#include "ITS_Bias.h"

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS ITS_BIAS
/*
Used to perform integrated tempering sampling (ITS) molecular dynamics
simulation on potential energy of the system or bias potential of CVs.

The ITS method generates a broad distribution as a function of the
potential energy \f$U\f$ of the system, which is achieved by using an
effective potential energy \f$U^{\text{eff}}\f$ :

\f[
U^{\text{eff}} = -\frac{1}{\beta_0} \log \sum_{k}^{N} n_k e^{- \beta_k U}
\f]

in which \f$\beta_0\f$ is the temperature of the system of interest,
\f$\beta_k\f$ is a series of temperatures that cover both low and high
temperatures, and \f$\{n_k\}\f$ are reweighting factors obtained through
an iterative procedure.

In PLUMED, all the bias are CV (collective variables) based. Therefore,
we change the form of ITS method to a CV-based bias potential \f$V(U)\f$:

\f[
V(U) = U^{\text{eff}} - U = -\frac{1}{\beta_0} \log \sum_{k}^{N} n_k e^{- \beta_k U} - U
\f]

In this form, the CV should be only chosen as energy, for example, the
potential energy of the system. In fact, the bias energy from the PLUMED
can be also used as the energy, in this cituation:

\f[
V(V') = -\frac{1}{\beta_0} \log \sum_{k}^{N} n_k e^{- \beta_0 (a_k V')} - V'
\f]

Therefore, the ITS method can be also consider about a broad distribution
as a summation of a series of rescaled bias potential. \f$a_k\f$ is the
rescale factor. In fact, \f$a_k\f$ and \f$\beta_k\f$ are equivalent because
\f$\beta_k = a_k * \beta_0\f$ then \f$a_k = T_0/T_k\f$. Therefore, you can choose the series of
temperatures or rescale factors to generate the broad distribution.

In order to obtain the weight factors \f$\{n_k\}\f$ at different temperature,
we use a iteration process during the ITS simulation. The principle of
setting the weight factors \f$\{n_k\}\f$ is let the distribution of potential
energy at different temperatures have the same contribution:

\f[
P_i=n_i \int_{r} e^{- \beta_i U} dr = n_j \int_{r} e^{- \beta_j U} dr = P_j
\f]
 
The system can be recovered to the thermodynamics at normal temperature
\f$\beta_0\f$ by multiplying a reweighting factor \f$c_0\f$ on the
probabilities of each observation:

\f[
c_0 = e^{\beta_0(U - U^{\text{eff}})} =
\left[\sum_{k}^{N} n_k e^{(\beta_0 - \beta_k)U}\right]^{-1}
\f]

\par Examples

\verbatim
energy: ENERGY

ITS_BIAS ...
  LABEL=its
  ARG=energy
  NREPLICA=100
  SIM_TEMP=300
  TEMP_MAX=370
  TEMP_MIN=270
  PACE=2000
  PESHIFT=1500
  FB_FILE=fb.data
  FB_STRIDE=100
  FBTRJ_FILE=fbtrj.data
  FBTRJ_STRIDE=20
... ITS_BIAS
\endverbatim

\verbatim
energy: ENERGY

DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 LABEL=restraint

ITS_BIAS ...
  LABEL=its
  ARG=energy
  BIAS=restraint
  NREPLICA=100
  SIM_TEMP=300
  TEMP_MAX=370
  TEMP_MIN=270
  PACE=2000
  PESHIFT=1500
  FB_FILE=fb.data
  FB_STRIDE=100
  FBTRJ_FILE=fbtrj.data
  FBTRJ_STRIDE=20
... ITS_BIAS
\endverbatim

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(ITS_Bias,"ITS_BIAS")

void ITS_Bias::registerKeywords(Keywords& keys)
{
	Bias::registerKeywords(keys);
	keys.addOutputComponent("rbias","default","the revised bias potential using rct");
	keys.addOutputComponent("energy","default","the instantaneous value of the potential energy of the system");
	keys.addOutputComponent("effective","DEBUG_FILE","the instantaneous value of the effective potential");
	keys.addOutputComponent("rct","DEBUG_FILE","the reweighting revise factor");
	keys.addOutputComponent("force","DEBUG_FILE","the instantaneous value of the bias force");
	keys.addOutputComponent("rwfb","DEBUG_FILE","the revised bias potential using rct");
	ActionWithValue::useCustomisableComponents(keys);
	keys.remove("ARG");
    keys.add("optional","ARG","the argument(CVs) here must be the potential energy. If no argument is setup, only bias will be integrated in this method."); 
	//~ keys.use("ARG");

	keys.add("compulsory","NREPLICA","the number of the replicas");
	keys.add("compulsory","PESHIFT","0.0","the shift value of potential energy");
	keys.add("compulsory","PACE","1000","the frequency for updating fb value");
	
	keys.add("optional","BIAS","the label of the bias to be used in the bias-ITS method");
	keys.add("optional","BIAS_RATIO","the ratio of the bias to be used in the bias-ITS method");
	keys.add("optional","SIM_TEMP","the temperature used in the simulation");
	keys.add("optional","TEMP_MIN","the lower bounds of the temperatures");
	keys.add("optional","TEMP_MAX","the upper bounds of the temperatures");
	keys.add("optional","REPLICA_RATIO_MIN","the minimal ratio of bias ratio when using IBS method");
	keys.add("optional","REPLICA_RATIO_MAX","the maximal ratio of bias ratio when using IBS method");
	keys.add("optional","TARGET_TEMP","the temperature that the system distribution at");
	
	keys.addFlag("EQUIVALENT_TEMPERATURE",false,"to simulate the system at target temperatue but keep using the original thermal bath");
	keys.addFlag("FB_FIXED",false,"to fix the fb value without update");
	keys.addFlag("MULTIPLE_WALKERS",false,"use multiple walkers");
	keys.addFlag("VARIATIONAL",false,"use variational approach to iterate the fb value");
	keys.addFlag("ONLY_FIRST_ORDER",false,"only to calculate the first order of Omega (gradient) in variational approach");
	keys.addFlag("TEMP_CONTRIBUTE",false,"use the contribution of each temperatue to calculate the derivatives instead of the target distribution");
	keys.addFlag("USE_FIXED_PESHIFT",false,"use the fixed PESHIFT value and do not update it automatically during the iteration of FB value");
	keys.addFlag("UNLINEAR_REPLICAS",false,"to setup the segments of temperature be propotional to the temperatues. If you setup the REPLICA_RATIO_MIN value, this term will be automatically opened.");
	keys.addFlag("DIRECT_AVERAGE",false,"to directly calculate the average of rbfb value in each step (only be used in traditional iteration process)");

	keys.add("optional","START_CYCLE","the start step for fb updating");
	keys.add("optional","FB_INIT","( default=0.0 ) the default value for fb initializing");
	keys.add("optional","RB_FAC1","( default=0.5 ) the ratio of the average value of rb");
	keys.add("optional","RB_FAC2","( default=0.0 ) the ratio of the old steps in rb updating");
	keys.add("optional","STEP_SIZE","( default=1.0 )the step size of fb iteration");

	keys.add("optional","FB_FILE","( default=fb.data ) a file to record the new fb values when they are update");
	keys.add("optional","FB_STRIDE","( default=1 ) the frequency to output the new fb values");
	keys.add("optional","FBTRAJ_STRIDE","( default=FB_STRIDE )the frequency to record the evoluation of fb values");
	keys.add("optional","FBTRAJ_FILE"," a file recording the evolution of fb values");
	//~ keys.add("optional","NORM_TRAJ"," a file recording the evolution of normalize factors");
	//~ keys.add("optional","ITER_TRAJ"," a file recording the evolution of the fb iteration factors");
	//~ keys.add("optional","DERIV_TRAJ"," a file recording the evolution of the derivation of fb factors");
	//~ keys.add("optional","PESHIFT_TRAJ"," a file recording the evolution of peshift");
	
	keys.add("optional","RW_TEMP","the temperatures used in the calcaulation of reweighting factors");
	keys.add("optional","RW_FILE","a file of the reweighting factors");
	keys.add("optional","RW_STRIDE","the frequency of reweighting factor output");
	keys.add("optional","RW_START","the start step of reweighting factor output");
	//~ keys.addFlag("REVISED_REWEIGHT",false,"calculate the c(t) value at reweighting factor output file");
	keys.add("optional","FB_READ_FILE","a file of reading fb values (include temperatures and peshift)");
	keys.add("optional","BIAS_FILE","a file to output the function of bias potential");
	keys.add("optional","RBFB_FILE","a file to output the evoluation of rbfb");
	keys.add("optional","BIAS_STRIDE","the frequency to output the bias potential");
	keys.add("optional","BIAS_MIN","the minimal value of coordinate of the bias potential function");
	keys.add("optional","BIAS_MAX","the maximal value of coordinate of the bias potential function");
	keys.add("optional","BIAS_BIN","the number of bins to record the bias potential function");
	keys.add("optional","DEBUG_FILE","a file of output debug information");
	keys.add("optional","POTDIS_FILE","the file that record the distribution of potential energy during the simulation");
	keys.add("optional","POTDIS_STEP","the freuency to record the distribution of potential energy");
	keys.add("optional","POTDIS_UPDATE","the freuency to output the distribution of potential energy");
	keys.add("optional","POTDIS_MIN","the minimal value of potential in the potential distribution file");
	keys.add("optional","POTDIS_MAX","the maximal value of potential in the potential distribution file");
	keys.add("optional","POTDIS_BIN","the number of bins that record the potential distribution");
	keys.add("optional","ENERGY_MIN","the minimal potential energy to calculate the 2nd order derivative");
	keys.add("optional","ENERGY_MAX","the maximal potential energy to calculate the 2nd order derivative");
	keys.add("optional","ENERGY_ACCURACY","the interval of potential energy to calculate the 2nd order derivative");

	keys.add("optional","ITERATE_LIMIT","to limit the iteration of fb value in variational in order to prevent the weight of higher temperatues become larger than the lower one");
}

ITS_Bias::~ITS_Bias()
{
	if(!equiv_temp&&!is_const)
	{
		ofb.close();
		ofbtrj.close();
		//~ onormtrj.close();
		//~ opstrj.close();
		//~ if(is_ves)
			//~ oderiv.close();
	}
	if(rw_output)
		ofw.close();
	if(is_debug)
		odebug.close();
	if(potdis_output)
		opotdis.close();
	if(rbfb_output)
		orbfb.close();
}

ITS_Bias::ITS_Bias(const ActionOptions& ao):
	PLUMED_BIAS_INIT(ao),update_start(0),rct(0),
	step(0),norm_step(0),mcycle(0),iter_limit(0),
	fb_init(0.0),fb_bias(0.0),rb_fac1(0.5),rb_fac2(0.0),step_size(1.0),
	is_const(false),rw_output(false),is_ves(false),
	read_norm(false),only_1st(false),bias_output(false),
	rbfb_output(false),is_debug(false),potdis_output(false),
	bias_linked(false),only_bias(false),
	is_set_temps(false),is_set_ratios(false),is_norm_rescale(false),
	read_fb(false),read_iter(false),fbtrj_output(false),
	rw_stride(1),rw_start(0),start_cycle(0),fb_stride(1),
	bias_stride(1),potdis_step(1),rctid(0),
	min_ener(1e38),pot_bin(1),dU(1),dvp2_complete(0)
{
	if(getNumberOfArguments()==0)
		only_bias=true;
	else if(getNumberOfArguments()>1)
		plumed_merror("this edition of ITS can only accept one CV");
		
	std::vector<std::string> bias_labels(0);
	parseVector("BIAS",bias_labels);
	nbiases = bias_labels.size();
	if(nbiases>0)
	{
		bias_linked=true;
		bias_pntrs_.resize(nbiases);
		for(unsigned int i=0; i<nbiases; i++)
		{
			bias_pntrs_[i]=plumed.getActionSet().selectWithLabel<Bias*>(bias_labels[i]);
			if(!bias_pntrs_[i]){plumed_merror("bias "+bias_labels[i]+" does not exist. NOTE: the Bias-ITS should always be defined AFTER all the other biases.");}
			bias_pntrs_[i]->linkExternalBias(this);
		}
		parseVector("BIAS_RATIO",bias_ratio);
		if(bias_ratio.size()==1)
		{
			double tbr=bias_ratio[0];
			bias_ratio.assign(nbiases,tbr);
		}
		else if(bias_ratio.size()!=nbiases)
			bias_ratio.resize(nbiases,1.0);
	}
	else if(only_bias)
		plumed_merror("the quantity of ARG and BIAS must be setup at last one");

	kB=plumed.getAtoms().getKBoltzmann();

	parseFlag("VARIATIONAL",is_ves);
	parseFlag("FB_FIXED",is_const);
	parseFlag("EQUIVALENT_TEMPERATURE",equiv_temp);
	parseFlag("MULTIPLE_WALKERS",use_mw);
	parseFlag("ONLY_FIRST_ORDER",only_1st);
	parseFlag("USE_FIXED_PESHIFT",use_fixed_peshift);
	parseFlag("UNLINEAR_REPLICAS",is_unlinear);
	parseFlag("DIRECT_AVERAGE",is_direct);
	parse("PESHIFT",peshift);
	parse("FB_READ_FILE",fb_input);
	
	if(only_bias&&equiv_temp)
		plumed_merror("EQUIVALENT_TEMPERATURE must be used with the potential energy as the argument");

	unsigned _nreplica(0);
	double _templ(-1),_temph(-1),_ratiol(-1),_ratioh(-1);
	parse("NREPLICA",_nreplica);
	if(_nreplica==0&&fb_input.size()==0)
		plumed_merror("NREPLICA must be set up or read from file.");
	parse("TEMP_MIN",_templ);
	parse("TEMP_MAX",_temph);
	if(_templ>0&&_temph>0)
		is_set_temps=true;
	parse("REPLICA_RATIO_MIN",_ratiol);
	parse("REPLICA_RATIO_MAX",_ratioh);
	
	if(_ratiol>0||_ratioh>0)
	{
		if(is_set_temps&&fb_input.size()==0)
			plumed_merror("the range of temperatures (TEMP_MIN and TEMP_MAX) and the range of replica ratio (REPLICA_RATIO_MIN and REPLICA_RATIO_MAX) cannot be setup simultaneously.");
		is_set_ratios=true;
		if(_ratiol<0)
		{
			if(is_unlinear)
				plumed_merror("if you want to use unliner temperatures, the minimal replica ratio (REPLICA_RATIO_MIN) must be setup!");
			_ratiol=0.0;
		}
		else
			is_unlinear=true;
		if(_ratioh<0) _ratioh=1.0;
	}
	if(!is_set_ratios&&!is_set_temps&&fb_input.size()==0)
		plumed_merror("the range of temperatures must be setup (TEMP_MIN/TEMP_MAX or REPLICA_RATIO_MIN/REPLICA_RATIO_MAX) or read from file (FB_READ_FILE).");
	
	parse("START_CYCLE",start_cycle);

	ener_min=peshift;
	parse("ENERGY_MIN",ener_min);
	ener_max=ener_min+1000;
	parse("ENERGY_MAX",ener_max);
	parse("ENERGY_ACCURACY",dU);

	double _kB=kB;
	double _peshift=peshift;
	
	sim_temp=-1;
	parse("SIM_TEMP",sim_temp);
	if(sim_temp>0)
		kBT=kB*sim_temp;
	else
	{
		kBT=plumed.getAtoms().getKbT();
		sim_temp=kBT/kB;
	}
	beta0=1.0/kBT;
	
	parse("START_CYCLE",start_cycle);
	
	unsigned read_count=0;
	
	if(fb_input.size()>0)
	{
		read_fb=true;
		read_count=read_fb_file(fb_input,_kB,_peshift);
		if(read_count==0)
			plumed_merror("do not read anything in the file "+fb_file);
		comm.Barrier();
		if(use_mw && comm.Get_rank()==0)
			multi_sim_comm.Barrier();
	}
	else
	{
		if(_nreplica<1)
			plumed_merror("the number of temperatures much be larger than 1");
		nreplica=_nreplica+1;

		if(is_set_ratios)
		{
			if(_ratiol<0)
				plumed_merror("REPLICA_RATIO_MIN must be large than 0");
			ratiol=_ratiol;

			if(_ratioh<0)
				plumed_merror("REPLICA_RATIO_MAX must be large than 0");
			ratioh=_ratioh;
			if(ratiol>=ratioh)
				plumed_merror("the value of RATIO_BIAS_MAX must be large than RATIO_BIAS_MIN");
				
			templ=sim_temp/ratioh;
			temph=sim_temp/ratiol;
			
			if(!is_unlinear)
			{
				double dr=(ratioh-ratiol)/_nreplica;
				for(unsigned i=0;i!=nreplica;++i)
				{
					double ratio_now=ratiol+i*dr;
					int_ratios.push_back(ratio_now);
					int_temps.push_back(sim_temp/ratio_now);
				}
			}
		}
		else
		{
			if(_templ<0)
				plumed_merror("TEMP_MIN must be large than 0");
			templ=_templ;

			if(_temph<0)
				plumed_merror("TEMP_MAX must be large than 0");
			temph=_temph;
			if(templ>=temph)
				plumed_merror("the value of TEMP_MAX must be large than TEMP_MIN");
				
			ratiol=sim_temp/temph;
			ratioh=sim_temp/templ;
			
			if(!is_unlinear)
			{
				double dt=(temph-templ)/_nreplica;
				for(unsigned i=0;i!=nreplica;++i)
				{
					double temp_now=templ+i*dt;
					int_temps.push_back(temp_now);
					int_ratios.push_back(sim_temp/temp_now);
				}
			}
		}
		if(is_unlinear)
		{
			double temp_ratio=exp(std::log(temph/templ)/_nreplica);
			for(unsigned i=0;i!=nreplica;++i)
			{
				double temp_now=templ*pow(temp_ratio,i);
				int_temps.push_back(temp_now);
				int_ratios.push_back(sim_temp/temp_now);
			}
		}
		
		mcycle=start_cycle;
		++start_cycle;
	}

	target_temp=sim_temp;
	parse("TARGET_TEMP",target_temp);
	if(equiv_temp)
	{
		kBT_target=kB*target_temp;
		beta_target=1.0/kBT_target;
		eff_factor=beta_target/beta0;
		bias_force=1-eff_factor;
	}
	//~ else if(sim_temp>temph||sim_temp<templ)
	//~ {
		//~ std::string sts,stl,sth;
		//~ Tools::convert(sim_temp,sts);
		//~ Tools::convert(templ,stl);
		//~ Tools::convert(temph,sth);
		//~ error("the value of SIM_TEMP("+sts+") must between TEMP_MIN("+
			//~ stl+") and TEMP_MAX("+sth+")");
	//~ }
	parse("FB_INIT",fb_init);
	parse("RB_FAC1",rb_fac1);
	parse("RB_FAC2",rb_fac2);
	parse("STEP_SIZE",step_size);
	
	if(fabs(step_size-1)>1.0e-6)
	{
		is_norm_rescale=true;
		fb_bias=std::log(step_size);
	}
	
	if(!read_norm)
		norml.assign(nreplica,0);

	gU.assign(nreplica,0);
	gf.assign(nreplica,0);
	bgf.assign(nreplica,0);
	rbfb.assign(nreplica,0);

	parse("PACE",update_step);
	if(update_step==0)
		plumed_merror("PACE cannot be 0");
		
	fb_file="fb.data";
	parse("FB_FILE",fb_file);
	parse("FB_STRIDE",fb_stride);

	parse("FBTRAJ_FILE",fb_trj);
	if(fb_trj.size()>0)
		fbtrj_output=true;
	fbtrj_stride=fb_stride;
	parse("FBTRAJ_STRIDE",fbtrj_stride);

	//~ parse("NORM_TRAJ",norm_trj);
	//~ parse("ITER_TRAJ",iter_trj);
	//~ if(((!is_ves)&&norm_trj.size()>0)||(is_ves&&iter_trj.size()>0))
		//~ norm_output=true;
	
	//~ deriv_trj="derivtrj.data";
	//~ parse("DERIV_TRAJ",deriv_trj);
	
	//~ parse("PESHIFT_TRAJ",peshift_trj);
	//~ if(peshift_trj.size()>0)
		//~ peshift_output=true;
	
	if(getRestart())
	{
		read_count=read_fb_file(fb_file,_kB,_peshift);
		if(read_count==0)
			plumed_merror("do not read anything in the file "+fb_file);
		comm.Barrier();
		if(use_mw && comm.Get_rank()==0)
			multi_sim_comm.Barrier();
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
			
		if(!read_fb&&!getRestart())
			fb.push_back(-exp(fb_init*(betak[0]-betak[i])));
	}
	
	rctid=find_rw_id(sim_temp,sim_dtl,sim_dth);
	fb0=find_rw_fb(rctid,sim_dtl,sim_dth);
	rct=-(fb0+std::log(double(nreplica)))/beta0;

	if(!equiv_temp&&!is_const)
	{
		setupOFile(fb_file,ofb,use_mw);
		ofb.addConstantField("ITERATE_METHOD");
		ofb.addConstantField("ITERATE_STEP");
		ofb.addConstantField("BOLTZMANN_CONSTANT");
		ofb.addConstantField("PESHIFT");
		ofb.addConstantField("COEFFICIENT_TYPE");
		ofb.addConstantField("NREPLICA");
		if(is_ves&&!only_1st)
		{
			ofb.addConstantField("ENERGY_MIN");
			ofb.addConstantField("ENERGY_MAX");
		}

		if(fbtrj_output)
			setupOFile(fb_trj,ofbtrj,use_mw);

		//~ if(is_ves)
		//~ {
			//~ if(norm_output)
				//~ setupOFile(iter_trj,onormtrj,use_mw);
			//~ 
			//~ setupOFile(deriv_trj,oderiv,use_mw);
		//~ }
		//~ else if(norm_output)
			//~ setupOFile(norm_trj,onormtrj,use_mw);
		
		//~ if(peshift_output)
			//~ setupOFile(peshift_trj,opstrj,use_mw);
	}
	
	parse("BIAS_FILE",bias_file);
	if(bias_file.size()>0)
	{
		bias_output=true;
		bias_max=0;
		parse("BIAS_MAX",bias_max);
		bias_min=0;
		parse("BIAS_MIN",bias_min);
		if(bias_max<=bias_min)
			plumed_merror("BIAS_MAX must be larger than bias_min");
		parse("BIAS_BIN",bias_bins);
		if(bias_bins<=1)
			plumed_merror("BIAS_BIN must be larger than 1");
		parse("BIAS_STRIDE",bias_stride);
		d_pot=(bias_max-bias_min)/(bias_bins-1);
	}
	
	parse("RBFB_FILE",rbfb_file);
	if(rbfb_file.size()>0)
	{
		rbfb_output=true;
		orbfb.link(*this);
		orbfb.open(rbfb_file);
	}
	
	parse("POTDIS_FILE",potdis_file);
	if(potdis_file.size()>0)
	{
		potdis_output=true;
		opotdis.link(*this);
		opotdis.open(potdis_file);
		parse("POTDIS_MIN",pot_min);
		parse("POTDIS_MAX",pot_max);
		parse("POTDIS_BIN",pot_bin);
		parse("POTDIS_STEP",potdis_step);
		potdis_update=update_step;
		parse("POTDIS_UPDATE",potdis_update);
		potdis_num=ceil((pot_max-pot_min)/pot_bin);
		pot_dis.assign(potdis_num,0);
	}

	parse("ITERATE_LIMIT",iter_limit);

	if(is_ves)
	{
		if(!is_const)
		{
			deriv_ps.assign(nreplica,1.0/nreplica);
			if(!only_1st)
			{
				if(ener_max<=ener_min)
				plumed_merror("TARGET_MAX must be larger than TARGET_MIN");

				ntarget=floor((ener_max-ener_min)/dU)+1;
				//~ if(!(read_mean&&read_dev))
					//~ plumed_merror("2nd order derivative calculation is open but there is no mean or deviation value in fb recording file");
				dvp2_complete=calc_deriv2();
				if(dvp2_complete<0.7)
				{
					std::string dpc;
					Tools::convert(dvp2_complete,dpc);
					plumed_merror("The compuatational comple of 2nd order derivative ("+dpc+") is too low!");
				}
			}
		}
		if(!read_iter)
		{
			double fb_zero=fb[0];
			for(unsigned i=0;i!=nreplica;++i)
			{
				fb[i]-=fb_zero;
				fb_iter.push_back(fb[i]);
			}
		}
	}

	parseVector("RW_TEMP",rw_temp);
	parse("RW_FILE",rw_file);
	parse("RW_STRIDE",rw_stride);
	parse("RW_START",rw_start);

	if(rw_file.size()==0&&rw_temp.size()>0)
		plumed_merror("RW_TEMP is setted but RW_FILE is not setted");

	if(rw_file.size()>0)
	{
		//~ parseFlag("REVISED_REWEIGHT",is_rw_rct);
		rw_output=true;
		if(rw_temp.size()>0)
		{
			rw_rctid.resize(rw_temp.size());
			rw_dth.resize(rw_temp.size());
			rw_dtl.resize(rw_temp.size());
			for(unsigned i=0;i!=rw_temp.size();++i)
			{
				if(rw_temp[i]<=temph&&rw_temp[i]>=templ)
					rw_beta.push_back(1.0/(kB*rw_temp[i]));
				else
					plumed_merror("the reweighting temperature must between TEMP_MIN and TEMP_MAX");
				rw_rctid[i]=find_rw_id(rw_temp[i],rw_dtl[i],rw_dth[i]);
			}
		}
		else
		{
			plumed_merror("RW_TEMP is not setup");
		}
		rw_factor.resize(rw_temp.size());

		setupOFile(rw_file,ofw,use_mw);
		for(unsigned i=0;i!=rw_temp.size();++i)
		{
			std::string tt;
			Tools::convert(i,tt);
			tt="rw_temperature_"+tt;
			ofw.addConstantField(tt).printField(tt,rw_temp[i]);
		}
	}

	parse("DEBUG_FILE",debug_file);
	if(debug_file.size()>0)
	{
		is_debug=true;
		odebug.link(*this);
		odebug.open(debug_file);
		
		odebug.addConstantField("ITERATE_STEP");
		odebug.addConstantField("NORM_STEP");
		odebug.addConstantField("STEP_SIZE");
		odebug.addConstantField("PESHIFT");
	}

	checkRead();

	addComponent("rbias"); componentIsNotPeriodic("rbias");
	valueRBias=getPntrToComponent("rbias");
	addComponent("energy"); componentIsNotPeriodic("energy");
	valueEnergy=getPntrToComponent("energy");
	
	if(is_debug)
	{
		addComponent("effective"); componentIsNotPeriodic("effective");
		valueEff=getPntrToComponent("effective");
		addComponent("rwfb"); componentIsNotPeriodic("rwfb");
		valueRwfb=getPntrToComponent("rwfb");
		addComponent("rct"); componentIsNotPeriodic("rct");
		valueRct=getPntrToComponent("rct");
		addComponent("force"); componentIsNotPeriodic("force");
		valueForce=getPntrToComponent("force");
		valueRwfb->set(fb0);
		valueRct->set(rct);
	}

	log.printf("  with simulation temperature: %f\n",sim_temp);
	log.printf("  with boltzmann constant: %f\n",kB);
	log.printf("  wiht beta (1/kT): %f\n",beta0);
	
	if(equiv_temp)
	{
		log.printf("  with equivalent target temperature: %f\n",target_temp);
		log.printf("  The system will run on target temperature (%fK) but keep the old thermal bath (%f).\n",target_temp,sim_temp);
	}
	else
	{
		if(only_bias)
			log.printf("  Only bias potential are setup as the CVs.\n");
		if(bias_linked)
		{
			if(nbiases==1)
				log.printf("  linked with bias \"%s\", and the ratio of bias is %f.\n",bias_labels[0].c_str(),bias_ratio[0]);
			else
				log.printf("  linked with %d biases:\n",nbiases);
			for(unsigned i=0;i!=nbiases;++i)
				log.printf("    Bias %d: %s with ratio %f\n",i+1,bias_labels[i].c_str(),bias_ratio[i]);
		}
		log.printf("  with basic parameters:\n");
		log.printf("    FB_INIT: %f\n",fb_init);
		log.printf("    FB_BIAS: %f\n",fb_bias);
		log.printf("    RB_FAC1: %f\n",rb_fac1);
		log.printf("    RB_FAC2: %f\n",rb_fac2);
		if(read_fb)
		{
			log.printf("  Reading in FB values from file: %s\n",fb_input.c_str());
			if(read_count==1)
				log.printf("  the reading iteration step is %d\n",int(mcycle));
			else
				log.printf("  with totally reading %d sets of coefficients, and the lastest one of iteration step %d will be used\n",int(read_count),int(mcycle));
		}
		else if(getRestart())
		{
			log.printf("  Restart running with reading in FB values from file: %s\n",fb_file.c_str());
			log.printf("  with totally reading %d sets of coefficients, and the lastest one of iteration step %d will be used\n",int(read_count),int(mcycle));
		}
		else if(is_unlinear)
			log.printf("  Unlinear replicas is open");

		if(read_norm)
		{
			log.printf("  with temperatue and FB value: (index temperature replica_ratio beta_k fb norml)\n");
			for(unsigned i=0;i!=nreplica;++i)
				log.printf("    %d\t%f\t%f\t%f\t%f\t%f\n",i,int_temps[i],int_ratios[i],betak[i],fb[i],norml[i]);
		}
		else
		{
			log.printf("  with temperatue and FB value: (index temperature replica_ratio beta_k fb)\n");
			for(unsigned i=0;i!=nreplica;++i)
				log.printf("    %d\t%f\t%f\t%f\t%f\n",i,int_temps[i],int_ratios[i],betak[i],fb[i]);
		}
		log.printf("    with temperatues of FB from %f to %f\n",templ,temph);
		log.printf("    with ratio of replica from %f to %f\n",ratiol,ratioh);
		
		log.printf("    with number of replica: %d\n",_nreplica);
		log.printf("    with PESHIFT value: %f\n",peshift);
		log.printf("    using the linear interpolation between the %d-th (%fK) and %d-th (%fK) FB value to calculate the reweight-revise factor\n",int(rctid),int_temps[rctid],int(rctid+1),int_temps[rctid+1]);
		if(fabs(_kB/kB-1)>1.0e-8)
			log.printf("    with Original PESHIFT value: %f (with boltzmann constant %f)\n",
				_peshift,_kB);
		if(use_mw)
		{
			log.printf("  Using multiple walkers");
			log.printf("   with number of walkers: %d\n",multi_sim_comm.Get_size());
			log.printf("   with walker number: %d\n",multi_sim_comm.Get_rank());
			log.printf("\n");
		}
		if(is_ves)
		{
			log.printf("  Using variational approach to iterate the fb value\n");
			if(!only_1st)
			{
				log.printf("    with calculation of Hessian is open\n");
				if(is_set_ratios)
				{
					for(unsigned i=0;i!=nreplica;++i)
						log.printf("    %d\t%f\t%f\t%f\n",i,int_ratios[i],energy_mean[i],energy_dev[i]);
				}
				else
				{
					for(unsigned i=0;i!=nreplica;++i)
						log.printf("    %d\t%f\t%f\t%f\n",i,int_temps[i],energy_mean[i],energy_dev[i]);
				}
				log.printf("    with energy integration range: from %f to %f\n",ener_min,ener_max);
				log.printf("    with integrative accuracy: %f\n",dU);
				log.printf("    with integrative complete: %f%\n",dvp2_complete*100);
				if(dvp2_complete<0.95)
					log.printf("    Warning! The computational complete of 2nd order derivative is low!");
			}
			log.printf("\n");
		}
		if(!is_const)
		{
			log.printf("  Using ITS update\n");
			if(is_direct)
				log.printf("    with directly calculating the avaerage of rbfb at each step\n");
			else
				log.printf("    with calculating the avaerage of rbfb at the end of fb updating\n");
			log.printf("    with frequence of FB value update: %d\n",
				update_step);
			log.printf("    writing FB output to file: %s\n",fb_file.c_str());
			log.printf("    writing FB trajectory to file: %s\n",fb_trj.c_str());
			if(rbfb_output)
				log.printf("    writing RBFB trajectory to file: %s\n",rbfb_file.c_str());
			log.printf("    writing normalized factors trajectory to file: %s\n",
				norm_trj.c_str());
			log.printf("    writing potential energy shift trajectory to file: %s\n",
				peshift_trj.c_str());
			log.printf("\n");
		}
		if(rw_output)
		{
			log.printf("  Writing reweighting factor to file: %s\n",rw_file.c_str());
			log.printf("    with reweighting factor at temperature:");
			for(unsigned i=0;i!=rw_temp.size();++i)
				log.printf(" %f",rw_temp[i]);
			log.printf("\n");
			log.printf("    with frequence of reweighting factor output: %d\n",rw_stride);
		}
		if(is_debug)
			log.printf("  Using debug mod with output file: %s\n",debug_file.c_str());

		log<<"Bibliography "<<
			plumed.cite("Gao, J. Chem. Phys. 128, 064105 (2008)");
		if(!is_ves)
			log<<plumed.cite("Gao, J. Chem. Phys. 128, 134111 (2008)");
		if(bias_linked)
			log<<plumed.cite("Yang, Zhang, Che, Yang, and Gao, J. Chem. Phys. 144, 094105 (2016)");
		log<<"\n";
	}
	if(bias_output)
	{
		log.printf("  Wirting output bias potential function\n");
		log.printf("    with output range from %f to %f\n",bias_min,bias_max);
		log.printf("    with output bins: %d\n",bias_bins);
		log.printf("    with Frequence of output: %d\n",bias_stride);
		log.printf("    wirting bias output to file: %s\n",bias_file.c_str());
		log.printf("\n");
	}
}

void ITS_Bias::calculate()
{
	if(only_bias)
		energy=0;
	else
		energy=getArgument(0);

	cv_energy=energy;
	
	if(bias_linked)
	{
		tot_bias=0;
		for(unsigned i=0;i!=nbiases;++i)
			tot_bias += bias_ratio[i] * bias_pntrs_[i]->getBias();
		cv_energy += tot_bias;
	}
	++step;

	// U = U_total + U_shift
	shift_energy=cv_energy+peshift;

	if(equiv_temp)
	{
		eff_energy=shift_energy*eff_factor;
		bias_energy=-1*shift_energy*bias_force;
	}
	else
	{
		for(unsigned i=0;i!=nreplica;++i)
		{
			// -\beta_k*U
			gU[i]=-betak[i]*shift_energy;
			// log[n_k*exp(-\beta_k*U)]
			gf[i]=gU[i]+fb[i];
			// log[\beta_k*n_k*exp(-\beta_k*U)]
			bgf[i]=gf[i]+std::log(betak[i]);
		}

		// log{\sum_k[n_k*exp(-\beta_k*U)]}
		gfsum=gf[0];
		// log{\sum_k[\beta_k*n_k*exp(-\beta_k*U)]}
		bgfsum=bgf[0];
		for(unsigned i=1;i!=nreplica;++i)
		{
			//~ exp_added(gUsum,gU[i]);
			exp_added(gfsum,gf[i]);
			exp_added(bgfsum,bgf[i]);
		}

		// U_EFF=-1/\beta_0*log{\sum_k[n_k*exp(-\beta_k*U)]}
		eff_energy=-gfsum/beta0;
		bias_energy=eff_energy-shift_energy;
		eff_factor=exp(bgfsum-gfsum)/beta0;
		bias_force=1.0-eff_factor;
	}

	setBias(bias_energy);
	valueRBias->set(bias_energy-rct);
	
	if(bias_linked)
	{
		for(unsigned i=0;i!=nbiases;++i)
			bias_pntrs_[i]->setExtraBiasRatio(-1.0*bias_force*bias_ratio[i]);
	}
	
	if(!only_bias)
		setOutputForce(0,bias_force);

	
	valueEnergy->set(shift_energy);
	
	if(is_debug)
	{
		valueForce->set(bias_force);
	}
	
	if(potdis_output)
	{
		if(step%potdis_step==0 && energy<pot_max && energy>pot_min)
		{
			unsigned id=floor((energy-pot_min)/pot_bin);
			++pot_dis[id];
		}
		if(step%potdis_update==0)
		{
			if(use_mw)
				multi_sim_comm.Sum(pot_dis);
			opotdis.rewind();
			opotdis.addConstantField("STEP").printField("STEP",int(step));
			opotdis.addConstantField("BOLTZMANN_CONSTANT").printField("BOLTZMANN_CONSTANT",kB);
			for(unsigned i=0;i!=potdis_num;++i)
			{
				double pot=pot_min+i*pot_bin;
				opotdis.printField("energy",pot);
				opotdis.printField("count",int(pot_dis[i]));
				opotdis.printField();
			}
			opotdis.flush();
		}
	}

	if(!equiv_temp)
	{
		if(!is_const)
		{
			if(is_ves||is_direct)
				update_rbfb_direct();
			else
				update_rbfb();
			++norm_step;
		}
		
		if(cv_energy<min_ener)
			min_ener=cv_energy;

		if(bias_output&&step==0)
			output_bias();

		if(rw_output && step%rw_stride==rw_start)
		{
			ofw.fmtField(" %f");
			ofw.printField("time",getTime());
			for(unsigned i=0;i!=rw_factor.size();++i)
			{
				std::string tt;
				Tools::convert(rw_temp[i],tt);
				rw_factor[i]=calc_bias(rw_beta[i]);
				//~ std::string tt1="BIAS_T"+tt;
				//~ ofw.printField(tt1,rw_factor[i]);
				double rwrct=-(find_rw_fb(rw_rctid[i],rw_dtl[i],rw_dth[i])+std::log(double(nreplica)))/rw_beta[i];
				std::string tt2="RBIAS_T"+tt;
				ofw.printField(tt2,rw_factor[i]-rwrct);
			}
			ofw.printField();
			ofw.flush();
		}

		if(is_debug)
		{
			if(!is_const)
			{
				energy_record.push_back(energy);
				gfsum_record.push_back(gfsum);
				bgfsum_record.push_back(bgfsum);
				effpot_record.push_back(eff_energy);
				bias_record.push_back(bias_energy);
				force_record.push_back(bias_force);
			}
			else
			{
				odebug.printField("step",int(step));
				odebug.printField("energy",energy);
				if(bias_linked)
					odebug.printField("tot_bias",tot_bias);
				odebug.printField("cv_energy",cv_energy);
				odebug.printField("peshift",peshift);
				odebug.printField("shift_energy",shift_energy);
				odebug.printField("eff_energy",eff_energy);
				odebug.printField("bias_energy",bias_energy);
				odebug.printField();
			}
		}

		if(!is_const && step%update_step==update_start && step>update_start)
		{
			if(is_debug)
			{
				for(unsigned i=0;i!=energy_record.size();++i)
				{
					int rs=step-energy_record.size()+i;
					odebug.printField("step",rs);
					odebug.printField("energy",energy_record[i]);
					odebug.printField("gfsum",gfsum_record[i]);
					odebug.printField("bgfsum",bgfsum_record[i]);
					odebug.printField("effective",effpot_record[i]);
					odebug.printField("bias",bias_record[i]);
					odebug.printField("force",force_record[i]);
					odebug.printField();
				}
				energy_record.resize(0);
				gfsum_record.resize(0);
				bgfsum_record.resize(0);
				effpot_record.resize(0);
				bias_record.resize(0);
				force_record.resize(0);
			}
			
			++mcycle;
			
			if(rbfb_output&&mcycle%fbtrj_stride==0)
			{
				orbfb.fmtField(" %f");
				orbfb.printField("step",int(mcycle));
				for(unsigned i=0;i!=nreplica;++i)
				{
					std::string id;
					Tools::convert(i,id);
					std::string rbfbid="RBFB"+id;
					orbfb.printField(rbfbid,rbfb[i]);
				}
				orbfb.printField();
				orbfb.flush();
			}

			//~ comm.Sum(rbfb);
			if(use_mw)
				mw_merge_rbfb();
			if(!use_fixed_peshift&&-1.0*min_ener>peshift)
				change_peshift(-1.0*min_ener);

			if(is_ves)
				fb_variational();
			else
				fb_iteration();
			
			fb0=find_rw_fb(rctid,sim_dtl,sim_dth);
			rct=-(fb0+std::log(double(nreplica)))/beta0;
			
			if(is_debug)
			{
				valueRwfb->set(fb0);
				valueRct->set(rct);
			}
				
			if(is_debug)
				odebug.flush();

			output_fb();
			for(unsigned i=0;i!=nreplica;++i)
				rbfb[i]=-1e38;

			norm_step=0;
		}
	}
}

// The iterate process of rbfb
// rbfb[k]=log[\sum_t(P_k)]
inline void ITS_Bias::update_rbfb()
{
	// the summation record the data of each the update steps (default=100)
	if(norm_step==0)
	{
		for(unsigned i=0;i!=nreplica;++i)
			rbfb[i]=gf[i];
	}
	else
	{
		for(unsigned i=0;i!=nreplica;++i)
			exp_added(rbfb[i],gf[i]);
	}
}

// The iterate process of rbfb
// rbfb[k]=log[\sum_t(p_k)]; p_k=P_k/[\sum_k(P_k)];
// The equivalence in variational iterate process:
// rbfb[i]=log[(\sum_t(\beta*(\partial V_bias(U;a)/\partial a_i)))_V(a)]
inline void ITS_Bias::update_rbfb_direct()
{
	//~ if(step%update_step==update_start)
	if(norm_step==0)
	{
		for(unsigned i=0;i!=nreplica;++i)
			rbfb[i]=gf[i]-gfsum;
	}
	else
	{
		for(unsigned i=0;i!=nreplica;++i)
			exp_added(rbfb[i],gf[i]-gfsum);
	}
}

void ITS_Bias::mw_merge_rbfb()
{
	multi_sim_comm.Sum(norm_step);
	
	unsigned nw=0;
	if(comm.Get_rank()==0)
		nw=multi_sim_comm.Get_size();
	comm.Bcast(nw,0);
	
	std::vector<double> all_min_ener(nw,0);
	std::vector<double> all_rbfb(nw*nreplica,0);
	if(comm.Get_rank()==0)
	{
		multi_sim_comm.Allgather(rbfb,all_rbfb);
		multi_sim_comm.Allgather(min_ener,all_min_ener);
	}
	comm.Bcast(all_rbfb,0);
	comm.Bcast(all_min_ener,0);

	min_ener=all_min_ener[0];
	for(unsigned j=0;j!=rbfb.size();++j)
		rbfb[j]=all_rbfb[j];
	for(unsigned i=1;i<nw;++i)
	{
		if(all_min_ener[i]<min_ener)
			min_ener=all_min_ener[i];
		for(unsigned j=0;j!=rbfb.size();++j)
			exp_added(rbfb[j],all_rbfb[i*rbfb.size()+j]);
	}
	if(is_ves||is_direct)
	{
		for(unsigned i=1;i!=nreplica;++i)
			rbfb[i]-=std::log(static_cast<double>(nw));
	}
}

// Y. Q. Gao, J. Chem. Phys. 128, 134111 (2008)
void ITS_Bias::fb_iteration()
{
	if(is_direct)
	{
		for(unsigned i=1;i!=nreplica;++i)
			rbfb[i]-=std::log(static_cast<double>(update_step));
	}
	else
	{
		double rbfbsum=rbfb[0];
		for(unsigned i=1;i!=nreplica;++i)
			exp_added(rbfbsum,rbfb[i]);
		for(unsigned i=0;i!=nreplica;++i)
			rbfb[i]-=rbfbsum;
	}
	
	// ratio[k]=log[m_k(t)]
	std::vector<double> ratio;
	if(mcycle==start_cycle&&!read_norm)
	{
		for(unsigned i=0;i!=nreplica-1;++i)
		{
			norml[i]=(rbfb[i]+rbfb[i+1])*rb_fac1;
			// m_k(0)=p_k(0)/p_{k+1}(0)
			ratio.push_back(rbfb[i+1]-rbfb[i]);
		}
	}
	else
	{
		for(unsigned i=0;i!=nreplica-1;++i)
		{
			// rb=log[f(p_k,p_{k+1};t)]
			// f(p_k,p_{k+1};t)=(p_k(t)*p_{k+1}(t))^c1+t*c2
			// (default f(p_k,p_{k+1};t)=\sqrt[p_k(t)*p_{k+1}(t)])
			double rb=(rbfb[i]+rbfb[i+1])*rb_fac1+(mcycle-1)*rb_fac2;
			// ratio_old=log[m_k(t-1)], m_k=n_k/n_{k+1}
			// (Notice that in the paper m_k=n_{k+1}/n_k)
			double ratio_old=fb[i]-fb[i+1];
			
			// normal=log[W_k(t)], normalold=log[W_k(t-1)]
			// W_k(t)=\sum_t[f(p_k,p_{k+1})]=W_k(t-1)+f(p_k,p_{k+1};t)
			// the summation record the data of ALL the update steps
			double normlold=norml[i];
			exp_added(norml[i],rb);
			
			// Rescaled normalization factor
			double rationorm=norml[i];
			// Wr_k(t)=W_k(t-1)+c_bias*f(p_k,p_{k+1};t)
			if(is_norm_rescale)
				rationorm=exp_add(normlold,rb+fb_bias);
				
			// m_k(t)=[c_bias*f(p_k,p_{k+1};t)*p_k(t)/p_{k+1}(t)+m_k(t-1)*W_k(t-1)]/Wr_k(t)
			ratio.push_back(exp_add(fb_bias+rb+rbfb[i+1]-rbfb[i],normlold+ratio_old)-rationorm);
		}
	}

	// partio[k]=log[1/n_k]
	std::vector<double> pratio(nreplica,0);
	// 1/n_k=n_1*\prod_j^k(m_j)
	for(unsigned i=0;i!=nreplica-1;++i)
		pratio[i+1]=pratio[i]+ratio[i];

	// As in the paper m_k=n_{k+1}/n_k
	// So here weights[k]=log[n_k]=log{1/pratio[k]}=-log{pratio[k]}
	
	//~ double rescalefb=std::log(double(nreplica))+fb[rctid]-betak[rctid]*peshift;
	
	bool is_nan=false;
	for(unsigned i=0;i!=nreplica;++i)
	{
		fb[i]=-pratio[i];
		if(fb[i]!=fb[i])
			is_nan=true;
	}
	if(is_debug||is_nan)
	{
		if(!is_debug)
		{
			odebug.link(*this);
			odebug.open("error.data");
		}
		
		odebug.printf("--- FB update ---\n");
		for(unsigned i=0;i!=nreplica;++i)
		{
			std::string id;
			Tools::convert(i,id);
			odebug.printField("index",id);
			if(is_set_ratios)
				odebug.printField("replica_ratio",int_ratios[i]);
			else
				odebug.printField("temperature",int_temps[i]);
			odebug.printField("betak",betak[i]);
			odebug.printField("rbfb",rbfb[i]);
			odebug.printField("fb",fb[i]);
			odebug.printField("norml",norml[i]);
			odebug.printField();
		}
		odebug.printf("--- FB update ---\n");
	}
	if(is_nan)
		plumed_merror("FB value become NaN. See debug file or \"error.data\" file for detailed information.");
}

void ITS_Bias::fb_variational()
{	
	// The average of 1st order drivative of bias potential in simulation
	std::vector<double> deriv_aver;
	// The Gradient of Omega
	std::vector<double> gradient;
	
	for(unsigned i=0;i!=nreplica;++i)
	{
		deriv_aver.push_back(exp(rbfb[i])/norm_step);
		gradient.push_back((deriv_aver[i]-deriv_ps[i])/beta0);
	}

	// fb_diff[i]=fb_iter[i]-fb[i]=a(i)-<a(i)>
	std::vector<double> fb_diff;
	// The Hessian of Omega
	std::vector<std::vector<double> > hessian(nreplica,std::vector<double>(nreplica,0));
	if(!only_1st)
	{
		std::vector<double>::const_iterator idv2=deriv2_ps.begin();
		for(unsigned i=0;i!=nreplica;++i)
		{
			fb_diff.push_back(fb_iter[i]-fb[i]);
			for(unsigned j=i;j!=nreplica;++j)
			{
				// <\partial V_bias(U;a)/\partial a_i>*<\partial V_bias(U;a)/\partial a_j>
				double d2a=deriv_aver[i]*deriv_aver[j];
				double d2p=*(idv2++);
				hessian[i][j]=(d2p-d2a)/beta0;
				if(i!=j)
					hessian[j][i]=hessian[i][j];
			}
		}
	}

	//~ double fb_zero;
	std::vector<double> o2as;
	//~ bool need_rescale=false;
	double fb_zero;
	bool is_nan=false;
	for(unsigned i=0;i!=nreplica;++i)
	{
		double fin_diff;
		if(only_1st)
			fin_diff=step_size*gradient[i];
		else
		{
			double o2a=0;
			for(unsigned j=0;j!=nreplica;++j)
				o2a+=hessian[i][j]*fb_diff[j];
			o2as.push_back(o2a);
			fin_diff=step_size*(gradient[i]+o2a);
		}

		if(i==0)
		{
			fb_zero=fin_diff;
			fin_diff=0;
		}
		else
		{
			fin_diff-=fb_zero;
			fb_iter[i]-=fin_diff;
			
			if(iter_limit>0)
			{
				unsigned itid;
				if(i>iter_limit)
					itid=i-iter_limit;
				else
					itid=0;
				
				if(fb_iter[i]>fb_iter[itid])
					fb_iter[i]=fb_iter[itid];
			}
			//~ else
			//~ {
				//~ double fb_i=fb[i]+(fb_iter[i]-fb[i])/(mcycle+1);
				//~ if(fb_i>=fb[i-1])
					//~ fb_iter[i]=fb_iter[i-1];
			//~ }

			fb[i]+=(fb_iter[i]-fb[i])/(mcycle+1);
		}
		if(fb[i]!=fb[i])
			is_nan=true;
	}

	//~ if(norm_output&&mcycle%fbtrj_stride==0)
	//~ {
		//~ oderiv.fmtField(" %f");
		//~ oderiv.printField("step",int(mcycle));
		//~ for(unsigned i=0;i!=nreplica;++i)
		//~ {
			//~ std::string id;
			//~ Tools::convert(i,id);
			//~ std::string fbid="FB"+id;
			//~ oderiv.printField(fbid,deriv_aver[i]);
		//~ }
		//~ oderiv.printField();
		//~ oderiv.flush();
	//~ }

	if(is_debug||is_nan)
	{
		if(!is_debug)
		{
			odebug.link(*this);
			odebug.open("error.data");
			odebug.addConstantField("ITERATE_STEP");
			odebug.addConstantField("NORM_STEP");
			odebug.addConstantField("STEP_SIZE");
			odebug.addConstantField("PESHIFT");
		}
		odebug.printf("--- FB Variational ---\n");
		odebug.printField("ITERATE_STEP",int(mcycle));
		odebug.printField("NORM_STEP",int(norm_step));
		odebug.printField("STEP_SIZE",step_size);
		odebug.printField("PESHIFT",peshift);
		for(unsigned i=0;i!=nreplica;++i)
		{
			odebug.printField("index",int(i));
			odebug.printField("rbfb",rbfb[i]);
			odebug.printField("deriv_aver",deriv_aver[i]);
			odebug.printField("deriv_ps",deriv_ps[i]);
			odebug.printField("gradient",gradient[i]);
			if(!only_1st)
				odebug.printField("o2a",o2as[i]);
			odebug.printField("fb",fb[i]);
			odebug.printField("fb_iter",fb_iter[i]);
			odebug.printField();
		}
		odebug.printf("--- Variational END ---\n");
	}
	if(is_nan)
		plumed_merror("FB value become NaN. See debug file or \"error.data\" file for detailed information.");
}

void ITS_Bias::output_fb()
{
	if(mcycle%fb_stride==0)
	{
		if(is_ves)
			ofb.printField("ITERATE_METHOD","Variational");
		else
			ofb.printField("ITERATE_METHOD","Traditional");
		ofb.printField("ITERATE_STEP",int(mcycle));
		ofb.printField("BOLTZMANN_CONSTANT",kB);
		ofb.printField("PESHIFT",peshift);
		ofb.printField("COEFFICIENT_TYPE",(is_set_temps?"TEMPERATURE":"RATIO"));
		ofb.printField("NREPLICA",int(nreplica));
		if(is_ves&&!only_1st)
		{
			ofb.printField("ENERGY_MIN",ener_min);
			ofb.printField("ENERGY_MAX",ener_max);
		}
		for(unsigned i=0;i!=nreplica;++i)
		{
			std::string id;
			Tools::convert(i,id);
			ofb.printField("index",id);
			if(is_set_ratios)
				ofb.printField("replica_ratio",int_ratios[i]);
			else
				ofb.printField("temperature",int_temps[i]);
			ofb.printField("fb_value",fb[i]);
			if(is_ves)
			{
				ofb.printField("fb_iteration",fb_iter[i]);
				if(!only_1st)
				{
					ofb.printField("energy_mean",energy_mean[i]);
					ofb.printField("energy_dev",energy_dev[i]);
				}
			}
			else
				ofb.printField("norm_value",norml[i]);
			ofb.printField();
		}
		ofb.printf("#!-----END-OF-FB-COEFFICIENTS-----\n\n");
		ofb.flush();
	}
	if(mcycle%fbtrj_stride==0)
	{
		if(fbtrj_output)
		{
			ofbtrj.fmtField(" %f");
			ofbtrj.printField("step",int(mcycle));
		}
		//~ if(norm_output)
		//~ {
			//~ onormtrj.fmtField(" %f");
			//~ onormtrj.printField("step",int(mcycle));
		//~ }
		if(fbtrj_output)
		{
			for(unsigned i=0;i!=nreplica;++i)
			{
				std::string id;
				Tools::convert(i,id);
				std::string fbid="FB"+id;
				
				if(fbtrj_output)
					ofbtrj.printField(fbid,fb[i]);

				//~ if(norm_output)
				//~ {
					//~ if(is_ves)
						//~ onormtrj.printField(fbid,fb_iter[i]);
					//~ else
						//~ onormtrj.printField(fbid,norml[i]);
				//~ }
			}
			ofbtrj.printField();
			ofbtrj.flush();
		}
		//~ if(norm_output)
		//~ {
			//~ onormtrj.printField();
			//~ onormtrj.flush();
		//~ }
		//~ if(peshift_output)
		//~ {
			//~ opstrj.fmtField(" %f");
			//~ opstrj.printField("step",int(mcycle));
			//~ opstrj.printField("PESHIFT",peshift);
			//~ opstrj.printField();
			//~ opstrj.flush();
		//~ }
	}
	if(bias_output&&mcycle%bias_stride==0)
		output_bias();
}

unsigned ITS_Bias::find_rw_id(double rwtemp,double& dtl,double& dth)
{
	if(rwtemp<templ)
		plumed_merror("the reweight temperature is lower than the minimal temperature of ITS");
	if(rwtemp>temph)
		plumed_merror("the reweight temperature is larger than the maximal temperature of ITS");
	
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
	
	return rwid;
}

double ITS_Bias::find_rw_fb(unsigned rwid,double dtl,double dth)
{
	double fbl=fb[rwid];
	double fbh=fb[rwid+1];
	double dfb=fbh-fbl;
	
	double rwfb;
	if(dtl<dth)
		rwfb=fbl+dfb*dtl/(dtl+dth);
	else
		rwfb=fbh-dfb*dth/(dtl+dth);
	return rwfb;
}

double ITS_Bias::find_rw_fb(double rwtemp)
{
	double dtl,dth;
	unsigned rwid=find_rw_id(rwtemp,dtl,dth);
	return find_rw_fb(rwid,dtl,dth);
}

void ITS_Bias::setupOFile(std::string& file_name, OFile& ofile, const bool multi_sim_single_files)
{
    ofile.link(*this);
    std::string fname=file_name;
    if(multi_sim_single_files)
    {
		unsigned int r=0;
		if(comm.Get_rank()==0)
			r=multi_sim_comm.Get_rank();
		comm.Bcast(r,0);
		if(r>0)
			fname="/dev/null";
		ofile.enforceSuffix("");
    }
    ofile.open(fname);
}

unsigned ITS_Bias::read_fb_file(const std::string& fname,double& _kB,double& _peshift)
{	
	IFile ifb;
    ifb.link(*this);
    if(use_mw)
		ifb.enforceSuffix("");
		
	if(!ifb.FileExist(fname))
		plumed_merror("Cannot find fb file " + fname );
		
    ifb.open(fname);
    
	ifb.allowIgnoredFields();
	
	std::string iter_method="Tradtional";
	
	unsigned _ntemp=nreplica;

	unsigned read_count=0;
	while(ifb)
	{
		if(ifb.FieldExist("ITERATE_METHOD"))
			ifb.scanField("ITERATE_METHOD",iter_method);

		if(ifb.FieldExist("ITERATE_STEP"))
		{
			int tmc=0;
			ifb.scanField("ITERATE_STEP",tmc);
			mcycle=tmc;
		}

		if(ifb.FieldExist("BOLTZMANN_CONSTANT"))
			ifb.scanField("BOLTZMANN_CONSTANT",_kB);

		if(ifb.FieldExist("PESHIFT"))
			ifb.scanField("PESHIFT",_peshift);
		
		std::string coe_type="TEMPERATURE";
		if(ifb.FieldExist("COEFFICIENT_TYPE"))
			ifb.scanField("COEFFICIENT_TYPE",coe_type);
		if(coe_type=="TEMPERATURE")
		{
			is_set_temps=true;
			is_set_ratios=false;
		}
		else if(coe_type=="RATIO")
		{
			is_set_temps=false;
			is_set_ratios=true;
		}
		else
			plumed_merror("unrecognized COEFFICIENT_TYPE "+coe_type);
			
		if(ifb.FieldExist("PESHIFT"))
			ifb.scanField("PESHIFT",_peshift);

		if(ifb.FieldExist("NREPLICA"))
		{
			int int_ntemp;
			ifb.scanField("NREPLICA",int_ntemp);
			_ntemp=int_ntemp;
		}
		
		if(is_set_ratios)
			int_ratios.resize(_ntemp);
		else
			int_temps.resize(_ntemp);
		fb.resize(_ntemp);

		if(ifb)
		{
			if(!ifb.FieldExist("fb_value"))
				plumed_merror("cannot found \"fb_value\" in file \""+fb_input+"\"");
			if(is_set_temps&&!ifb.FieldExist("temperature"))
				plumed_merror("cannot found \"temperature\" in file \""+fb_input+"\"");
			if(is_set_ratios&&!ifb.FieldExist("replica_ratio"))
				plumed_merror("cannot found \"replica_ratio\" in file \""+fb_input+"\"");
		}
		else
			break;

		bool read_mean=false;
		bool read_dev=false;

		if(!is_const)
		{
			if(is_ves)
			{
				if(ifb.FieldExist("fb_iteration"))
					read_iter=true;
				if(!only_1st)
				{
					if(ifb.FieldExist("energy_mean"))
						read_mean=true;
					if(ifb.FieldExist("energy_dev"))
						read_dev=true;
					if(ifb.FieldExist("ENERGY_MIN"))
						ifb.scanField("ENERGY_MIN",ener_min);
					if(ifb.FieldExist("ENERGY_MAX"))
						ifb.scanField("ENERGY_MAX",ener_max);
				}
			}
			else
			{
				if(ifb.FieldExist("norm_value"))
					read_norm=true;
			}
		}
		
		if(is_ves)
		{
			if(read_iter)
				fb_iter.resize(_ntemp);
			if(read_mean)
				energy_mean.resize(_ntemp);
			if(read_dev)
				energy_dev.resize(_ntemp);
		}
		else if(read_norm)
			norml.resize(_ntemp);
		
		double tmpfb;
		
		for(unsigned i=0;i!=_ntemp;++i)
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
			
			if(is_ves)
			{
				if(read_iter)
				{
					double tmpif;
					ifb.scanField("fb_iteration",tmpif);
					fb_iter[i]=tmpif;
				}
				if(read_mean)
				{
					double tmpmean;
					ifb.scanField("energy_mean",tmpmean);
					energy_mean[i]=tmpmean;
				}
				if(read_dev)
				{
					double tmpdev;
					ifb.scanField("energy_dev",tmpdev);
					energy_dev[i]=tmpdev;
				}
			}
			else if(read_norm)
			{
				double tmpnb;
				ifb.scanField("norm_value",tmpnb);
				norml[i]=tmpnb;
			}
			fb[i]=tmpfb;
			ifb.scanField();
		}
		++read_count;
	}
	if(iter_method=="Traditional")
	{
		if(is_ves)
			plumed_merror("The iterate method is variational but the fb input file is a traditional input file");
	}
	else if(iter_method=="Variational")
	{
		if(!is_ves)
			plumed_merror("The iterate method is traditional but the fb input file is a variational input file");
	}
	else
		plumed_merror("unknown iterate method: \""+iter_method+"\"");
	
	peshift=_peshift;
	
	if(is_set_ratios)
	{
		nreplica=_ntemp;
		ratiol=ratioh=int_ratios[0];
		for(unsigned i=1;i!=int_ratios.size();++i)
		{
			if(int_ratios[i]>ratioh)
				ratioh=int_ratios[i];
			if(int_ratios[i]<ratiol)
				ratiol=int_ratios[i];
		}
		temph=sim_temp/ratiol;
		templ=sim_temp/ratioh;
	}
	else
	{
		nreplica=_ntemp;
		templ=temph=int_temps[0];
		for(unsigned i=1;i!=int_temps.size();++i)
		{
			if(int_temps[i]>temph)
				temph=int_temps[i];
			if(int_temps[i]<templ)
				templ=int_temps[i];
		}
		ratioh=sim_temp/templ;
		ratiol=sim_temp/temph;
	}
	if(fabs(_kB/kB-1)>1.0e-8)
	{
		double rescale=kB/_kB;
		peshift*=rescale;
		if(!only_1st)
		{
			ener_max*=rescale;
			ener_min*=rescale;
			for(unsigned i=0;i!=nreplica;++i)
			{
				energy_mean[i]*=rescale;
				energy_dev[i]*=rescale;
			}
		}
	}
	ifb.close();
	return read_count;
}

void ITS_Bias::change_peshift(double new_shift)
{
	fb_rescale(new_shift-peshift);
	if(is_ves)
		iter_rescale(new_shift-peshift);
	peshift=new_shift;
}

void ITS_Bias::output_bias()
{
	std::string smcycle;
	Tools::convert(int(mcycle),smcycle);
	std::string file=bias_file+smcycle+".data";
	
	obias.link(*this);
	obias.open(file);
	
	if(step==0)
	{
		obias.addConstantField("ITERATE_STEP");
		obias.addConstantField("BOLTZMANN_CONSTANT");
		obias.addConstantField("PESHIFT");
	}

	obias.printField("ITERATE_STEP",int(mcycle));
	obias.printField("BOLTZMANN_CONSTANT",kB);
	obias.printField("PESHIFT",peshift);

	double sfb=fb[0];
	for(unsigned i=1;i!=nreplica;++i)
		exp_added(sfb,fb[i]);

	for(unsigned i=0;i!=bias_bins;++i)
	{
		double energy=bias_min+d_pot*i;
		double senergy=energy+peshift;

		double _gfsum=-betak[0]*senergy+fb[0];
		double _bgfsum=_gfsum+std::log(betak[0]);
		for(unsigned j=1;j!=nreplica;++j)
		{
			double tmp=-betak[j]*senergy+fb[j];
			exp_added(_gfsum,tmp);
			exp_added(_bgfsum,tmp+std::log(betak[j]));
		}

		double eff_pot=-_gfsum/beta0+sfb;
		double bias_ener=eff_pot-senergy;
		double fscale=exp(_bgfsum-_gfsum)/beta0;
		obias.printField("original_potential",energy);
		obias.printField("bias_potential",bias_ener);
		obias.printField("effective_potential",eff_pot);
		obias.printField("force_scale",fscale);
		obias.printField();
	}
	obias.flush();
	obias.close();
}

double ITS_Bias::calc_deriv2()
{
	deriv2_ps.assign(nreplica*(nreplica+1)/2,0);
	double complete=0;
	for(unsigned id_U=0;id_U!=ntarget;++id_U)
	{
		double ener=ener_min+id_U*dU;
		std::vector<double> prob;
		double psum=0;
		for(unsigned i=0;i!=nreplica;++i)
		{
			double pp=exp(-pow((ener-energy_mean[i])/energy_dev[i],2)/2)/energy_dev[i];
			prob.push_back(pp);
			psum+=pp;
		}
		complete+=psum;
		std::vector<double>::iterator idv2=deriv2_ps.begin();
		for(unsigned i=0;i!=nreplica;++i)
		{
			for(unsigned j=i;j!=nreplica;++j)
			{
				double vv=(prob[i]*prob[j]/psum);
				*(idv2++)+=vv;
			}
		}
	}
	double normsum=dU/(nreplica*sqrt(2*pi));
	complete*=normsum;

	for(unsigned i=0;i!=deriv2_ps.size();++i)
		deriv2_ps[i]*=normsum;
	return complete;
}

std::vector<double> ITS_Bias::omega2_alpha(const std::vector<std::vector<double> >& hessian,const std::vector<double>& dalpha)
{
	std::vector<double> result;
	for(unsigned i=0;i!=hessian.size();++i)
	{
		double tmp=0;
		for(unsigned j=0;j!=dalpha.size();++j)
			tmp+=hessian[i][j]*dalpha[j];
		result.push_back(tmp);
	}
	return result;
}

std::vector<double> ITS_Bias::omega2_alpha(const std::vector<double>& hessian,const std::vector<double>& dalpha)
{
	std::vector<double> result;
	for(unsigned i=0;i!=dalpha.size();++i)
		result.push_back(hessian[i]*dalpha[i]);
	return result;
}

}
}
