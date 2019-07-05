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
	keys.addOutputComponent("rct","default","the reweighting revise factor");
	keys.addOutputComponent("force","default","the instantaneous value of the bias force");
	keys.addOutputComponent("potential","default","the instantaneous value of the potential energy of the system");
	keys.addOutputComponent("effective","DEBUG_FILE","the instantaneous value of the effective potential");
	//~ keys.addOutputComponent("rbias_T","RW_TEMP","the revised bias potential at different temperatures");
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
	keys.add("optional","RATIO_MIN","the minimal ratio of temperature");
	keys.add("optional","RATIO_MAX","the maximal ratio of temperature");
	keys.add("optional","TARGET_TEMP","the temperature that the system distribution at");

	keys.addFlag("EQUIVALENT_TEMPERATURE",false,"to simulate the system at target temperatue but keep using the original thermal bath");
	keys.addFlag("FB_FIXED",false,"to fix the fb value without update");
	keys.addFlag("MULTIPLE_WALKERS",false,"use multiple walkers");
	keys.addFlag("VARIATIONAL",false,"use variational approach to iterate the fb value");
	keys.addFlag("ONLY_FIRST_ORDER",false,"only to calculate the first order of Omega (gradient) in variational approach");
	keys.addFlag("TEMP_CONTRIBUTE",false,"use the contribution of each temperatue to calculate the derivatives instead of the target distribution");
	keys.addFlag("PESHIFT_AUTO_ADJUST",false,"automatically adjust the PESHIFT value during the fb iteration");
	keys.addFlag("UNLINEAR_REPLICAS",false,"to setup the segments of temperature be propotional to the temperatues. If you setup the REPLICA_RATIO_MIN value, this term will be automatically opened.");
	//~ keys.addFlag("DIRECT_AVERAGE",false,"to directly calculate the average of rbfb value in each step (only be used in traditional iteration process)");
	keys.addFlag("NOT_USE_BIAS_RCT",false,"do not use the c(t) of bias to modify the bias energy when using bias as the input CVs");
	keys.addFlag("REWEIGHT_PARTITION",false,"use reweighting at different temperatures to calculate partition functions");
	keys.addFlag("ABSOLUTE_RATIO",false,"use the absolute ratio of the partition functions to calculate the fb values");

	keys.add("optional","START_CYCLE","the start step for fb updating");
	keys.add("optional","FB_INIT","( default=0.0 ) the default value for fb initializing");
	keys.add("optional","RB_FAC1","( default=0.5 ) the ratio of the average value of rb");
	keys.add("optional","RB_FAC2","( default=0.0 ) the ratio of the old steps in rb updating");
	keys.add("optional","FB_BIAS","( default=0.0 ) the step size of fb iteration");
	keys.add("optional","TEMP_RATIO_ENERGY","( default=0 ) the energy to adjust the ratio of each temperatures");
	keys.add("optional","ANNEAL_RATE","( default=0 ) the rate of annealing per step. The value is the KL divergence of the new distribution with the old one of potential energy");
	keys.add("optional","RCT_FILE","the file to output the c(t)");
	keys.add("optional","RCT_STRIDE","the frequency to output the c(t)");
	keys.add("optional","ANNEAL_STRIDE","the frequncy (how many update cycles) to anneal");
	keys.add("optional","EMA_DAYS","( default=0 ) the periods for smoothing factor of the exponential moving average to iteration the fb value");
	//~ keys.add("optional","PRE_ITERATION_CYCLE","( default=100 ) the periods for smoothing factor of the exponential moving average to iteration the fb value");

	keys.add("optional","FB_FILE","( default=fb.data ) a file to record the new fb values when they are update");
	keys.add("optional","FB_STRIDE","( default=1 ) the frequency to output the new fb values");
	keys.add("optional","FBTRAJ_STRIDE","( default=FB_STRIDE )the frequency to record the evoluation of fb values");
	keys.add("optional","FBTRAJ_FILE"," a file recording the evolution of fb values");
	//~ keys.add("optional","NORM_TRAJ"," a file recording the evolution of normalize factors");
	//~ keys.add("optional","ITER_TRAJ"," a file recording the evolution of the fb iteration factors");
	//~ keys.add("optional","DERIV_TRAJ"," a file recording the evolution of the derivation of fb factors");
	//~ keys.add("optional","PESHIFT_TRAJ"," a file recording the evolution of peshift");

	keys.add("optional","RW_TEMP","the temperatures used in the calcaulation of reweighting factors");
	//~ keys.addFlag("REVISED_REWEIGHT",false,"calculate the c(t) value at reweighting factor output file");
	keys.add("optional","FB_READ_FILE","a file of reading fb values (include temperatures and peshift)");
	keys.add("optional","BIAS_FILE","a file to output the function of bias potential");
	keys.add("optional","RBFB_FILE","a file to output the evoluation of rbfb");
	keys.add("optional","ZKE_FILE","a file to output the evoluation of the effective partition function Z_k(E)");
	keys.add("optional","ZKU_FILE","a file to output the evoluation of the partition function Z_k(U) of the system");
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
	keys.add("optional","PEGAUSS_FILE","the file that estimate the distribution of potential energy using gaussian distribution");
	keys.add("optional","PEGAUSS_STEP","( default=100 ) the freuency of update cycles to calculate the gaussian estimated distribution of potential energy");
	keys.add("optional","PEGAUSS_BINS","( default=1000 ) the number of bins to output the gaussian estimated distribution of potential energy");
	keys.add("optional","PEGAUSS_SHIFT","( default=0 ) the constant value to shift the potential energy at the integration of gaussian");
	keys.add("optional","GAUSS_INTEGRATE_RANGE","( default=5 ) the range (the times of SD of the Gaussian) to integrate the gaussian esitimated distribution of potential energy");
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
		if(rct_output)
			orct.close();
		//~ onormtrj.close();
		//~ opstrj.close();
	}
	if(is_debug)
		odebug.close();
	if(potdis_output)
		opotdis.close();
	if(rbfb_output)
		orbfb.close();
	if(zke_output)
		ozke.close();
	if(zku_output)
		ozku.close();
	if(do_calc_pegauss)
		opegauss.close();
}

ITS_Bias::ITS_Bias(const ActionOptions& ao):
	PLUMED_BIAS_INIT(ao),update_start(0),rct(0),step(0),norm_step(0),
	mcycle(0),iter_limit(0),ema_days(0),
	fb_init(0.0),fb_bias(0.0),rb_fac1(0.5),rb_fac2(0.0),temp_ratio_energy(0),
	is_const(false),rw_output(false),read_norm(false),only_1st(false),
	bias_output(false),rbfb_output(false),zke_output(false),
	zku_output(false),is_debug(false),potdis_output(false),
	bias_linked(false),only_bias(false),is_read_ratio(false),
	is_set_temps(false),is_set_ratios(false),is_norm_rescale(false),
	read_fb(false),read_iter(false),fbtrj_output(false),
	rct_output(false),read_zbu(false),rw_partition(false),
	abs_ratio(false),do_calc_pegauss(false),
	start_cycle(0),fb_stride(1),bias_stride(1),potdis_step(1),
	min_ener(1e38),pot_bin(1),dU(1),avg_shift(0),gauss_limit(5)
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

	parseFlag("FB_FIXED",is_const);
	parseFlag("EQUIVALENT_TEMPERATURE",equiv_temp);
	parseFlag("MULTIPLE_WALKERS",use_mw);
	parseFlag("ONLY_FIRST_ORDER",only_1st);
	parseFlag("PESHIFT_AUTO_ADJUST",auto_peshift);
	parseFlag("UNLINEAR_REPLICAS",is_unlinear);
	parseFlag("NOT_USE_BIAS_RCT",no_bias_rct);
	parseFlag("ABSOLUTE_RATIO",abs_ratio);
	parseFlag("REWEIGHT_PARTITION",rw_partition);
	parse("PESHIFT",peshift);
	parse("FB_READ_FILE",fb_input);
	fb_file="fb.data";
	parse("FB_FILE",fb_file);
	parse("EMA_DAYS",ema_days);
	//~ parse("PRE_ITERATION_CYCLE",pre_iter_steps);
	
	anneal_rate=0;
	parse("ANNEAL_RATE",anneal_rate);
	if(anneal_rate>1e-38)
	{
		do_calc_pegauss=true;
		do_anneal=true;
	}
	
	double wema=0;
	if(ema_days>0)
	{
		use_ema=true;
		wema=2.0/ema_days;
		wemat=std::log(wema);
		wemay=std::log(1.0-wema);
	}
	else
		use_ema=false;

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
	parse("RATIO_MIN",_ratiol);
	parse("RATIO_MAX",_ratioh);
	//~ parse("NEIGHBOR_RATIO",neighbor_ratio);
	parse("TEMP_RATIO_ENERGY",temp_ratio_energy);

	parse("PACE",update_step);
	if(update_step==0)
		plumed_merror("PACE cannot be 0");

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
	
	unsigned read_count=0;
	
	if(fb_input.size()>0)
	{
		read_fb=true;
		read_count=read_fb_file(fb_input,_kB,_peshift);
		if(read_count==0)
			plumed_merror("do not read anything in the file "+fb_input);
		comm.Barrier();
		if(use_mw && comm.Get_rank()==0)
			multi_sim_comm.Barrier();
	}
	else if(getRestart())
	{
		read_count=read_fb_file(fb_file,_kB,_peshift);
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
				
			templ=sim_temp*ratiol;
			temph=sim_temp*ratioh;
			
			if(!is_unlinear)
			{
				double dr=(ratioh-ratiol)/_nreplica;
				for(unsigned i=0;i!=nreplica;++i)
				{
					double ratio_now=ratiol+i*dr;
					int_ratios.push_back(ratio_now);
					int_temps.push_back(sim_temp*ratio_now);
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
				
			ratiol=templ/sim_temp;
			ratioh=temph/sim_temp;
			
			if(!is_unlinear)
			{
				double dt=(temph-templ)/_nreplica;
				for(unsigned i=0;i!=nreplica;++i)
				{
					double temp_now=templ+i*dt;
					int_temps.push_back(temp_now);
					int_ratios.push_back(temp_now/sim_temp);
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
				int_ratios.push_back(temp_now/sim_temp);
			}
		}
		mcycle=start_cycle;
		++start_cycle;
	}
	for(unsigned i=0;i!=nreplica;++i)
		betak.push_back(1.0/(kB*int_temps[i]));

	logN=std::log(double(nreplica));
	ratio_norm=logN;
	fb_ratio0=-logN;
	// the target distribution of P'_k=Z_k/\sum_i Z_i, default value is 1/N
	fb_ratios.assign(nreplica,-logN);
	// the target ratio of the two neighor distribution R_k=P'_{k+1}/P'_{k}, default value is 1:1
	rbfb_ratios.assign(nreplica,0);

	fb_rct.assign(nreplica,0);

	if(fabs(temp_ratio_energy)>1.0e-15)
		update_ratios(temp_ratio_energy);
	else
		do_anneal=false;

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
	parse("FB_BIAS",fb_bias);
	
	if(fabs(fb_bias>1e-8))
		is_norm_rescale=true;
	
	if(!read_norm)
	{
		norml.assign(nreplica,1e-38);
		zbe.assign(nreplica,1e-38);
	}
	
	if(!read_zbu)
		zbu.assign(nreplica,1e-38);

	gU.assign(nreplica,0);
	gE.assign(nreplica,0);
	gf.assign(nreplica,0);
	rbfb.assign(nreplica,1e-38);
	//~ rbfb2.assign(nreplica,1e-38);
	rbzb.assign(nreplica,1e-38);
	bgf.assign(nreplica,0);
	peshift_ratio.assign(nreplica,0);

	parse("FB_STRIDE",fb_stride);

	parse("FBTRAJ_FILE",fb_trj);
	if(fb_trj.size()>0)
		fbtrj_output=true;
	fbtrj_stride=fb_stride;
	parse("FBTRAJ_STRIDE",fbtrj_stride);
	
	parse("RCT_FILE",rct_file);
	if(rct_file.size()>0)
		rct_output=true;
	rct_stride=fb_stride;
	parse("RCT_STRIDE",rct_stride);
	
	// s(k)=exp(\beta_k*E_shift)/exp(\beta_{k+1}*E_shift)
	// s(k)=exp[(\beta_k-\beta_{k+1})*E_shift]
	set_peshift_ratio();
	if(!read_fb&&!getRestart())
	{
		for(unsigned i=0;i!=nreplica;++i)
		{
			//~ fb.push_back(-exp(fb_init*(betak[0]-betak[i])));
			fb.push_back(fb_init*(int_temps[i]-int_temps[0]));
		}
	}
	
	for(unsigned i=0;i!=nreplica;++i)
		fb_rct[i]=calc_rct(betak[i],fb[i],fb_ratios[i]);

	if(!equiv_temp&&!is_const)
	{
		setupOFile(fb_file,ofb,use_mw);
		ofb.addConstantField("ITERATE_METHOD");
		ofb.addConstantField("BOLTZMANN_CONSTANT");
		ofb.addConstantField("ITERATE_STEP");
		ofb.addConstantField("PACE");
		ofb.addConstantField("PESHIFT");
		ofb.addConstantField("COEFFICIENT_TYPE");
		ofb.addConstantField("NREPLICA");
		ofb.addConstantField("TEMP_RATIO_ENERGY");

		if(fbtrj_output)
			setupOFile(fb_trj,ofbtrj,use_mw);
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
		setupOFile(rbfb_file,orbfb,use_mw);
	}

	parse("ZKE_FILE",zke_file);
	if(zke_file.size()>0)
	{
		zke_output=true;
		ozke.link(*this);
		setupOFile(zke_file,ozke,use_mw);
	}

	parse("ZKU_FILE",zku_file);
	if(zku_file.size()>0)
	{
		zku_output=true;
		ozku.link(*this);
		setupOFile(zku_file,ozku,use_mw);
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
	
	parse("PEGAUSS_FILE",pegauss_file);
	output_pegauss=false;
	pegauss_update=1;
	dE_num=1000;
	if(pegauss_file.size()>0)
	{
		do_calc_pegauss=true;
		output_pegauss=true;
		avg_step=0;
		parse("PEGAUSS_STEP",pegauss_update);
		parse("PEGAUSS_BINS",dE_num);
		parse("PEGAUSS_SHIFT",avg_shift);
		parse("GAUSS_INTEGRATE_RANGE",gauss_limit);
		
		int_num=dE_num+1;
		
		setupOFile(pegauss_file,opegauss,use_mw);
		opegauss.addConstantField("ITERATION");
		opegauss.addConstantField("BOLTZMANN_CONSTANT");
		opegauss.addConstantField("TEMP_RATIO_ENERGY");
		
		peavg.assign(nreplica,1e-38);
		peavg2.assign(nreplica,1e-38);
		penorm.assign(nreplica,1e-38);
		
		gauss_mean.assign(nreplica,0);
		gauss_sd.assign(nreplica,0);
		gauss_sd2.assign(nreplica,0);
		gauss_weight0.assign(nreplica,0);
		gauss_weight1.assign(nreplica,0);
		gmin_id.resize(nreplica);
		gmax_id.resize(nreplica);
		
		energy_range.assign(int_num,0);
		eff_p0.assign(int_num,0);
		eff_p1.assign(int_num,0);
	}

	parse("ITERATE_LIMIT",iter_limit);
	
	rctid=find_rw_id(sim_temp);
	fb0=find_rw_fb(rctid);
	rct=calc_rct(beta0,fb0,fb_ratio0);
	
	addComponent("rbias"); componentIsNotPeriodic("rbias");
	valueRBias=getPntrToComponent("rbias");
	//~ addComponent("energy"); componentIsNotPeriodic("energy");
	//~ valueEnergy=getPntrToComponent("energy");
	//~ addComponent("rct"); componentIsNotPeriodic("rct");
	//~ valueRct=getPntrToComponent("rct");
	//~ valueRct->set(rct);
	setRctComponent("rct");
	setRct(rct);
	addComponent("force"); componentIsNotPeriodic("force");
	valueForce=getPntrToComponent("force");
	addComponent("Potential"); componentIsNotPeriodic("Potential");
	valuePot=getPntrToComponent("Potential");
	addComponent("effective"); componentIsNotPeriodic("effective");
	valueEff=getPntrToComponent("effective");

	parseVector("RW_TEMP",rw_temp);

	std::vector<double> rw_fb(rw_temp.size());
	if(rw_temp.size()>0)
	{
		rw_output=true;

		rw_rct.resize(rw_temp.size());
		rw_rctid.resize(rw_temp.size());
		rw_fb_ratios.resize(rw_temp.size());
		rw_factor.resize(rw_temp.size());
		valueRwbias.resize(rw_temp.size());
		for(unsigned i=0;i!=rw_temp.size();++i)
		{
			if(rw_temp[i]<=temph&&rw_temp[i]>=templ)
				rw_beta.push_back(1.0/(kB*rw_temp[i]));
			else
				plumed_merror("the reweighting temperature must between TEMP_MIN and TEMP_MAX");
			rw_rctid[i]=find_rw_id(rw_temp[i]);
			
			std::string tt;
			Tools::convert(rw_temp[i],tt);
			std::string rtt="rbias_T"+tt;
			addComponent(rtt); componentIsNotPeriodic(rtt);
			valueRwbias[i]=getPntrToComponent(rtt);
		}
		
		update_rw_rct(true);
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
	
	if(is_debug)
	{
		addComponent("rwfb"); componentIsNotPeriodic("rwfb");
		valueRwfb=getPntrToComponent("rwfb");
		valueRwfb->set(fb0);
	}

	checkRead();
	
	log<<"  with bibliography ";
	log<<plumed.cite("Gao, J. Chem. Phys. 128, 064105 (2008)");
	log<<plumed.cite("Yang, Niu and Parrinello, J. Phys. Chem. Lett. 9, 6426 (2018)");
	log<<"\n";

	log.printf("  with simulation temperature: %f\n",sim_temp);
	log.printf("  with boltzmann constant: %f\n",kB);
	log.printf("  with beta (1/kT): %f\n",beta0);
	log.printf("  with energy to adjust the ratio of temperatures: %f\n",temp_ratio_energy);
	log.printf("  with fb ratio at simulation temperature (%f): %f\n",sim_temp, fb_ratio0);
	
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
			if(no_bias_rct)
				log.printf("    without using the c(t) of the bias\n");
			else
				log.printf("    using the c(t) of the bias if necessary.\n");
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
			if(abs_ratio)
			{
				log.printf("  with temperatue and FB value: (index T_k gamma_k beta_k P'_k fb Z_b(E)\n");
				for(unsigned i=0;i!=nreplica;++i)
					log.printf("    %d\t%f\t%f\t%f\t%f\t%f\t%f\n",i,int_temps[i],int_ratios[i],betak[i],fb_ratios[i],fb[i],zbe[i]);
			}
			else
			{
				log.printf("  with temperatue and FB value: (index T_k gamma_k beta_k P'_k fb norml)\n");
				for(unsigned i=0;i!=nreplica;++i)
					log.printf("    %d\t%f\t%f\t%f\t%f\t%f\t%f\n",i,int_temps[i],int_ratios[i],betak[i],fb_ratios[i],fb[i],norml[i]);
			}
		}
		else
		{
			log.printf("  with temperatue and FB value: (index[k] T_k gamma_k beta_k P'_k fb_k)\n");
			for(unsigned i=0;i!=nreplica;++i)
				log.printf("    %d\t%f\t%f\t%f\t%f\t%f\n",i,int_temps[i],int_ratios[i],betak[i],fb_ratios[i],fb[i]);
		}
		log.printf("    with temperatues of FB from %f to %f\n",templ,temph);
		log.printf("    with ratio of replica from %f to %f\n",ratiol,ratioh);
		
		log.printf("    with number of replica: %d\n",_nreplica);
		log.printf("    with PESHIFT value: %f\n",peshift);
		log.printf("    using the linear interpolation between the fb values at %f-th temperature\n",rctid);
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

		if(!is_const)
		{
			log.printf("  Updating fb values\n");
			if(abs_ratio)
				log.printf("    using the absolute ratio of the partition function at each temperature\n");
			else
			{
				log.printf("    using the relative ratio between the partition functions between two neighbor temperatures.");
				log<<plumed.cite("Gao, J. Chem. Phys. 128, 134111 (2008)");
				log.printf("\n");
			}
				
			if(rw_partition)
				log.printf("    with integration of the direct avaerage to calculate the partition functions (Z_k(E)=\\sum_t(n_k*exp(-beta_k*E(t))))\n");
			else
				log.printf("    with integration of the reweight avaerage to calculate the partition functions (Z_k(E)=\\sum_t[exp(beta_k*{V_k[E(t)]-c_k(t)})])\n");
			log.printf("    with calculating the avaerage of rbfb at the end of fb updating\n");
			log.printf("    with frequence of FB value update: %d\n",update_step);
			if(use_ema)
				log.printf("    using exponential moving average with %d days (smoothing factor is %f)\n",int(ema_days),wema);
			else
				log.printf("    using running average.\n");
			log.printf("    writing FB output to file: %s\n",fb_file.c_str());
			log.printf("    writing FB trajectory to file: %s\n",fb_trj.c_str());
			if(rbfb_output)
				log.printf("    writing RBFB trajectory to file: %s\n",rbfb_file.c_str());
			log.printf("\n");
		}
		if(rw_output)
		{
			log.printf("  Output reweighting factor at temperature:\n");
			for(unsigned i=0;i!=rw_temp.size();++i)
				log.printf("    %d\t%fK(%f) fit at %f temperature with target ratio %f\n",int(i),rw_temp[i],rw_beta[i],rw_rctid[i],rw_fb_ratios[i]);
		}
		
		if(do_calc_pegauss)
		{
			log.printf("  Estimating the distribution of potential energy as Gaussian distribuition\n");
			log.printf("    with updating for each %d cycles\n",int(pegauss_update));
			log.printf("    with the range of integration as %f times of standard deviation\n",gauss_limit);
			log.printf("    with integration bins: %d\n",int(int_num));
			if(output_pegauss)
				log.printf("    with distribution output to file: %s\n",pegauss_file.c_str());
			if(do_anneal&&!is_const)
			{
				log.printf("  Using annealing with rate %f\n",anneal_rate);
			}
		}
		
		if(is_debug)
			log.printf("  Using debug mod with output file: %s\n",debug_file.c_str());
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

	shift_pot=energy+peshift;

	cv_energy=energy;

	if(bias_linked)
	{
		tot_bias=0;
		for(unsigned i=0;i!=nbiases;++i)
		{
			tot_bias += bias_ratio[i] * bias_pntrs_[i]->getBias();
			if((!no_bias_rct)&&bias_pntrs_[i]->isSetRct())
				tot_bias -= bias_pntrs_[i]->getRct();
		}
		cv_energy += tot_bias;
	}
	++step;

	// U = U_total + E_shift
	input_energy=cv_energy+peshift;

	if(equiv_temp)
	{
		eff_energy=input_energy*eff_factor;
		bias_energy=-1*input_energy*bias_force;
	}
	else
	{
		for(unsigned i=0;i!=nreplica;++i)
		{
			// -\beta_k*U_pot
			gU[i]=-betak[i]*energy;
			// -\beta_k*U
			gE[i]=-betak[i]*input_energy;
			// log[n_k*exp(-\beta_k*U)]
			gf[i]=gE[i]+fb[i];
			// log[\beta_k*n_k*exp(-\beta_k*U)]
			bgf[i]=gf[i]+std::log(betak[i]);
		}

		// log{\sum_k[n_k*exp(-\beta_k*U)]}
		gfsum=gf[0];
		// log{\sum_k[\beta_k*n_k*exp(-\beta_k*U)]}
		bgfsum=bgf[0];
		for(unsigned i=1;i!=nreplica;++i)
		{
			//~ exp_added(gUsum,gE[i]);
			exp_added(gfsum,gf[i]);
			exp_added(bgfsum,bgf[i]);
		}

		// U_EFF=-1/\beta_0*log{\sum_k[n_k*exp(-\beta_k*U)]}
		eff_energy=-gfsum/beta0;
		bias_energy=eff_energy-input_energy;
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

	//~ valueEnergy->set(input_energy);
	valuePot->set(shift_pot);
	valueEff->set(eff_energy);
	valueForce->set(eff_factor);

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
			update_rbfb();
			++norm_step;
		}
		if(do_calc_pegauss)
			update_peavg();
		
		if(cv_energy<min_ener)
			min_ener=cv_energy;

		if(bias_output&&step==0)
			output_bias();

		if(rw_output)
		{
			for(unsigned i=0;i!=rw_factor.size();++i)
			{
				rw_factor[i]=calc_bias(rw_beta[i]);
				valueRwbias[i]->set(rw_factor[i]-rw_rct[i]);
			}
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
				odebug.printField("input_energy",input_energy);
				odebug.printField("eff_energy",eff_energy);
				odebug.printField("bias_energy",bias_energy);
				odebug.printField();
			}
		}

		if(step%update_step==update_start && step>=update_step)
		{
			++mcycle;
			
			if(use_mw)
				mw_merge_data();
			if(auto_peshift&&-1.0*min_ener>peshift)
				change_peshift(-1.0*min_ener);
				
			zbu_iteration();
			if(!is_const)
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

				if(abs_ratio)
					fb_partition();
				else
					fb_iteration();
				
				fb0=find_rw_fb(rctid);
				rct=calc_rct(beta0,fb0,fb_ratio0);
				setRct(rct);
				
				for(unsigned i=0;i!=nreplica;++i)
					fb_rct[i]=calc_rct(betak[i],fb[i],fb_ratios[i]);
					
				if(rw_output)
					update_rw_rct(false);

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
				
				if(zke_output&&mcycle%fbtrj_stride==0)
				{
					ozke.fmtField(" %f");
					ozke.printField("step",int(mcycle));
					for(unsigned i=0;i!=nreplica;++i)
					{
						std::string id;
						Tools::convert(i,id);
						std::string zkeid="zke"+id;
						ozke.printField(zkeid,zbe[i]);
					}
					ozke.printField();
					ozke.flush();
				}
				
				if(zku_output&&mcycle%fbtrj_stride==0)
				{
					ozku.fmtField(" %f");
					ozku.printField("step",int(mcycle));
					for(unsigned i=0;i!=nreplica;++i)
					{
						std::string id;
						Tools::convert(i,id);
						std::string zkuid="zku"+id;
						ozku.printField(zkuid,zbu[i]);
					}
					ozku.printField();
					ozku.flush();
				}
				
				if(is_debug)
				{
					valueRwfb->set(fb0);
				}
					
				if(is_debug)
					odebug.flush();
			}
			
			if(do_calc_pegauss&&mcycle%pegauss_update==0)
			{
				calc_pe_vag();
				update_gaussian_grid();
				
				if(!is_const&&do_anneal)
				{
					annealing_fb();

					fb0=find_rw_fb(rctid);
					rct=calc_rct(beta0,fb0,fb_ratio0);
					setRct(rct);
					
					for(unsigned i=0;i!=nreplica;++i)
						fb_rct[i]=calc_rct(betak[i],fb[i],fb_ratios[i]);
						
					if(rw_output)
						update_rw_rct(true);
				}
				
				if(output_pegauss)
				{
					opegauss.printField("ITERATION",int(mcycle));
					opegauss.printField("BOLTZMANN_CONSTANT",kB);
					opegauss.printField("TEMP_RATIO_ENERGY",temp_ratio_energy);
					for(unsigned i=0;i!=int_num;++i)
					{
						opegauss.printField("energy",energy_range[i]);
						opegauss.printField("probability",eff_p0[i]);
						opegauss.printField();
					}
					opegauss.flush();
				}
				
				avg_step=0;
			}
			
			if(!is_const)
			{
				output_fb();
			}
			
			for(unsigned i=0;i!=nreplica;++i)
			{
				rbfb[i]=-1e38;
				rbzb[i]=-1e38;
			}

			norm_step=0;
		}
	}
}


inline void ITS_Bias::update_rbfb()
{	
	// the summation record the data of each the update steps (default=100)
	if(rw_partition)
	{
		if(norm_step==0)
		{
			for(unsigned i=0;i!=nreplica;++i)
			{
				rbfb[i]=gE[i]-gfsum-betak[i]*fb_rct[i];
				rbzb[i]=gU[i]-gfsum-betak[i]*fb_rct[i];
			}
		}
		else
		{
			for(unsigned i=0;i!=nreplica;++i)
			{
				exp_added(rbfb[i],gE[i]-gfsum-betak[i]*fb_rct[i]);
				exp_added(rbzb[i],gU[i]-gfsum-betak[i]*fb_rct[i]);
			}
		}
	}
	else
	{
		// The iterate process of rbfb
		// rbfb[k]=log[\sum_t(Z_k)]
		if(norm_step==0)
		{
			for(unsigned i=0;i!=nreplica;++i)
			{
				rbfb[i]=gf[i];
				//~ rbfb2[i]=gE[i]-gfsum-betak[i]*fb_rct[i];
				rbzb[i]=gU[i]-gfsum-betak[i]*fb_rct[i];
			}
		}
		else
		{
			// The iterate process of rbfb
			// rbfb[k]=log[\sum_t(P_k)]; P_k=Z_k/[\sum_k(Z_k)];
			// The equivalence in variational iterate process:
			// rbfb[i]=log[(\sum_t(\beta*(\partial V_bias(U;a)/\partial a_i)))_V(a)]
			for(unsigned i=0;i!=nreplica;++i)
			{
				exp_added(rbfb[i],gf[i]);
				//~ exp_added(rbfb2[i],gE[i]-gfsum-betak[i]*fb_rct[i]);
				exp_added(rbzb[i],gU[i]-gfsum-betak[i]*fb_rct[i]);
			}
		}
	}
}

inline void ITS_Bias::update_peavg()
{
	double input_avg=input_energy+avg_shift;
	if(input_avg<=0)
	{
		double shiftdiff=-input_avg+1.0;
		double logdiff=std::log(shiftdiff);
		input_avg=1.0;
		// E'_i = E_i + c
		avg_shift+=shiftdiff;
		if(avg_step!=0)
		{
			for(unsigned i=0;i!=nreplica;++i)
			{
				// <(E_i + c)^2> = <E_i^2> + 2*c*<E_i> + c^2
				double avg2diff=exp_add(std::log(2.0)+logdiff+peavg[i],2*logdiff+penorm[i]);
				exp_added(peavg2[i],avg2diff);

				// <E_i + c> = <E_i> + c
				exp_added(peavg[i],logdiff+penorm[i]);
			}
		}
	}
	double logpe=std::log(input_avg);
	double logpe2=2*logpe;

	if(avg_step==0)
	{
		for(unsigned i=0;i!=nreplica;++i)
		{
			double weight=gE[i]-gfsum-betak[i]*fb_rct[i];
			peavg[i]=logpe+weight;
			peavg2[i]=logpe2+weight;
			penorm[i]=weight;
		}
	}
	else
	{
		for(unsigned i=0;i!=nreplica;++i)
		{
			double weight=gE[i]-gfsum-betak[i]*fb_rct[i];
			exp_added(peavg[i],logpe+weight);
			exp_added(peavg2[i],logpe2+weight);
			exp_added(penorm[i],weight);
		}
	}
	++avg_step;
}

inline void ITS_Bias::update_rw_rct(bool use_new_ratios)
{
	if(use_new_ratios)
	{
		for(unsigned i=0;i!=rw_temp.size();++i)
		{
			double rw_fb=find_rw_fb(rw_rctid[i]);
			//~ rw_fb_ratios[i]=rw_rctid[i]*rbfb_ratio-ratio_norm;
			rw_fb_ratios[i]=-temp_ratio_energy*rw_beta[i]-ratio_norm;
			rw_rct[i]=calc_rct(rw_beta[i],rw_fb,rw_fb_ratios[i]);
		}
	}
	else
	{
		for(unsigned i=0;i!=rw_temp.size();++i)
		{
			double rw_fb=find_rw_fb(rw_rctid[i]);
			rw_rct[i]=calc_rct(rw_beta[i],rw_fb,rw_fb_ratios[i]);
		}
	}
}

inline void ITS_Bias::update_ratios(double new_energ)
{
	double ratio1=-new_energ*betak[0];
	ratio_norm=ratio1;
	fb_ratios[0]=ratio1;
	for(unsigned i=1;i!=nreplica;++i)
	{
		double ratio_value=-new_energ*betak[i];
		fb_ratios[i]=ratio_value;
		exp_added(ratio_norm,ratio_value);
		rbfb_ratios[i-1]=fb_ratios[i]-fb_ratios[i-1];
	}
	fb_ratio0=-new_energ*beta0-ratio_norm;
	
	for(unsigned i=0;i!=nreplica;++i)
		fb_ratios[i]-=ratio_norm;
}

void ITS_Bias::calc_pe_vag()
{
	for(unsigned i=0;i!=nreplica;++i)
	{
		peavg[i]-=penorm[i];
		peavg2[i]-=penorm[i];
		
		double mean=exp(peavg[i]);
		double mean2=exp(peavg2[i]);
		gauss_mean[i]=mean-avg_shift-peshift;
		gauss_sd2[i]=mean2-mean*mean;
		gauss_sd[i]=sqrt(gauss_sd2[i]);
		gauss_weight0[i]=exp(fb_ratios[i])/(gauss_sd[i]*std::log(2.0*pi));
	}
}

void ITS_Bias::update_gaussian_grid()
{
	double eff_min=gauss_mean.front()-gauss_limit*gauss_sd.front();
	double eff_max=gauss_mean.back()+gauss_limit*gauss_sd.back();
	
	int_dE=(eff_max-eff_min)/dE_num;
	kl_rate=anneal_rate/int_dE;
	
	for(unsigned i=0;i!=energy_range.size();++i)
		energy_range[i]=eff_min+int_dE*i;
	for(unsigned i=0;i!=nreplica;++i)
	{
		double gmin=gauss_mean[i]-gauss_limit*gauss_sd[i];
		double gmax=gauss_mean[i]+gauss_limit*gauss_sd[i];
		unsigned id0=0;
		if(gmin>eff_min)
			id0=floor((gmin-eff_min)/int_dE);
		gmin_id[i]=id0;
		unsigned id1=ceil((gmax-eff_min)/int_dE)+1;
		if(id1>int_num+1)
			id1=int_num+1;
		gmax_id[i]=id1;
	}
	
	for(unsigned i=0;i!=eff_p0.size();++i)
		eff_p0[i]=0;

	gauss_values.resize(0);
	for(unsigned i=0;i!=nreplica;++i)
	{
		std::vector<double> tmp;
		for(unsigned j=gmin_id[i];j!=gmax_id[i];++j)
		{
			double e_avg=energy_range[j]-gauss_mean[i];
			double gauss=exp(-0.5*e_avg*e_avg/gauss_sd2[i]);
			tmp.push_back(gauss);
			eff_p0[j]+=gauss_weight0[i]*gauss;
		}
		gauss_values.push_back(tmp);
	}
}

bool ITS_Bias::anneal_judge(double new_energy)
{
	update_ratios(new_energy);
	
	for(unsigned i=0;i!=eff_p1.size();++i)
		eff_p1[i]=0;

	for(unsigned i=0;i!=nreplica;++i)
	{
		 gauss_weight1[i]=exp(fb_ratios[i])/(gauss_sd[i]*std::log(2.0*pi));
		 unsigned id=0;
		 for(unsigned j=gmin_id[i];j!=gmax_id[i];++j)
		 {
			 double gauss=gauss_values[i][id++];
			 eff_p1[j]+=gauss_weight1[i]*gauss;
		 }
	}
	
	return kl_divergence(eff_p0,eff_p1);
}

inline bool ITS_Bias::kl_divergence(const std::vector<double>& p0,const std::vector<double>& p1)
{
	double kl0=0;
	double kl1=0;
	for(unsigned i=0;i!=energy_range.size();++i)
	{
		kl0+=-p0[i]*std::log(p1[i]/p0[i]);
		kl1+=-p1[i]*std::log(p0[i]/p1[i]);
		if(kl0>kl_rate&&kl1>kl_rate)
			return true;
	}
	
	return false;
}

void ITS_Bias::annealing_fb()
{
	std::vector<double> old_ratios(fb_ratios);
	
	double new_energy=0;
	do_anneal=false;
	while(anneal_judge(new_energy))
	{
		do_anneal=true;
		new_energy+=(temp_ratio_energy-new_energy)/2;
	}
	
	temp_ratio_energy=new_energy;
	for(unsigned i=0;i!=nreplica;++i)
		fb[i]+=fb_ratios[i]-old_ratios[i];
	double fb0=fb[0];
	for(unsigned i=0;i!=nreplica;++i)
		fb[i]-=fb0;
}

void ITS_Bias::mw_merge_data()
{
	multi_sim_comm.Sum(norm_step);
	
	unsigned nw=0;
	if(comm.Get_rank()==0)
		nw=multi_sim_comm.Get_size();
	comm.Bcast(nw,0);
	
	std::vector<double> all_min_ener(nw,0);
	std::vector<double> all_rbfb(nw*nreplica,0);
	std::vector<double> all_rbzb(nw*nreplica,0);
	
	std::vector<double> all_peavg(nw*nreplica,0);
	std::vector<double> all_peavg2(nw*nreplica,0);
	std::vector<double> all_penorm(nw*nreplica,0);
	
	if(comm.Get_rank()==0)
	{
		multi_sim_comm.Allgather(rbfb,all_rbfb);
		multi_sim_comm.Allgather(rbzb,all_rbzb);
		multi_sim_comm.Allgather(min_ener,all_min_ener);
		if(do_calc_pegauss)
		{
			multi_sim_comm.Allgather(peavg,all_peavg);
			multi_sim_comm.Allgather(peavg2,all_peavg2);
			multi_sim_comm.Allgather(penorm,all_penorm);
		}
	}
	comm.Bcast(all_rbfb,0);
	comm.Bcast(all_rbzb,0);
	comm.Bcast(all_min_ener,0);
	if(do_calc_pegauss)
	{
		comm.Bcast(all_peavg,0);
		comm.Bcast(all_peavg2,0);
		comm.Bcast(all_penorm,0);
	}

	min_ener=all_min_ener[0];
	for(unsigned j=0;j!=nreplica;++j)
	{
		rbfb[j]=all_rbfb[j];
		rbzb[j]=all_rbzb[j];
		if(do_calc_pegauss)
		{
			peavg[j]=all_peavg[j];
			peavg2[j]=all_peavg2[j];
			penorm[j]=all_penorm[j];
		}
	}

	for(unsigned i=1;i<nw;++i)
	{
		if(all_min_ener[i]<min_ener)
			min_ener=all_min_ener[i];
		for(unsigned j=0;j!=nreplica;++j)
		{
			exp_added(rbfb[j],all_rbfb[i*nreplica+j]);
			exp_added(rbzb[j],all_rbzb[i*nreplica+j]);
			if(do_calc_pegauss)
			{
				exp_added(peavg[j],all_peavg[i*nreplica+j]);
				exp_added(peavg2[j],all_peavg2[i*nreplica+j]);
				exp_added(penorm[j],all_penorm[i*nreplica+j]);
			}
		}
	}

	for(unsigned i=0;i!=nreplica;++i)
	{
		rbfb[i]-=std::log(static_cast<double>(nw));
		rbzb[i]-=std::log(static_cast<double>(nw));
		if(do_calc_pegauss)
		{
			peavg[i]-=std::log(static_cast<double>(nw));
			peavg2[i]-=std::log(static_cast<double>(nw));
			penorm[i]-=std::log(static_cast<double>(nw));
		}
	}
}

// Y. Q. Gao, J. Chem. Phys. 128, 134111 (2008)
void ITS_Bias::fb_iteration()
{
	for(unsigned i=0;i!=nreplica;++i)
		rbfb[i]-=std::log(static_cast<double>(update_step));
	
	double rbfbsum=rbfb[0];
	for(unsigned i=1;i!=nreplica;++i)
		exp_added(rbfbsum,rbfb[i]);
	for(unsigned i=0;i!=nreplica;++i)
		rbfb[i]-=rbfbsum;

	// ratio[k]=log[m_k(t)]
	std::vector<double> ratio;
	std::vector<double> old_fb;
	if(is_debug)
		old_fb=fb;

	if(mcycle==start_cycle)
	{
		for(unsigned i=0;i!=nreplica-1;++i)
		{
			norml[i]=(rbfb[i]+rbfb[i+1])*rb_fac1;
			// m_k(0)=n_k(0)/n_{k+1}(0)
			double ratio_old=fb[i]-fb[i+1];
			// m_k(1)={m'}_k(1)=m_k(0)*P_k(0)/P_{k+1}(0)/R_k(0)
			ratio.push_back(ratio_old+rbfb[i+1]-rbfb[i]-rbfb_ratios[i]);
		}
	}
	else
	{
		for(unsigned i=0;i!=nreplica-1;++i)
		{
			// rb=log[f(P_k,P_{k+1};t)]
			// f(P_k,P_{k+1};t)=(P_k(t)*P_{k+1}(t))^c1+t*c2
			// (default f(P_k,P_{k+1};t)=\sqrt[P_k(t)*P_{k+1}(t)])
			double rb=(rbfb[i]+rbfb[i+1])*rb_fac1+(mcycle-1)*rb_fac2;
			// ratio_old=log[m_k(t-1)], m_k=n_k/n_{k+1}
			// (Notice that m_k in the paper equal to n_{k+1}/n_k)
			double ratio_old=fb[i]-fb[i+1];
			// ratio_new=log[{m'}_k(t)]=log[m_k(t-1)*P_{k+1}/P_{k}//R_k], if P_{k+1}/P_{k}=R_k, m_k(t)=m_k(t-1)
			double ratio_new=ratio_old+rbfb[i+1]-rbfb[i]-rbfb_ratios[i];
			
			if(use_ema&&mcycle>=ema_days)
				norml[i]+=wemay;

			// normal=log[W_k(t)], normalold=log[W_k(t-1)]
			// W_k(t)=\sum_t[f(P_k,P_{k+1})]=W_k(t-1)+f(P_k,P_{k+1};t)
			// the summation record the data of ALL the update steps
			double normlold=norml[i];
			exp_added(norml[i],rb);

			// Rescaled normalization factor
			double rationorm=norml[i];
			// the weight of new ratio
			// w_k(t)=c_bias*f(P_k,P_{k+1};t)
			double weight=fb_bias+rb;
			// Wr_k(t)=W_k(t-1)+w_k(t)
			if(is_norm_rescale)
				rationorm=exp_add(normlold,weight);

			// m_k(t)=[m_k(t-1)*c_bias*f(P_k,P_{k+1};t)*P_k(t)/P_{k+1}(t)+m_k(t-1)*W_k(t-1)]/Wr_k(t)
			// m_k(t)=[{m'}_k(t)*w_k(t)+m_k(t-1)*W_k(t-1)]/Wr_k(t)
			ratio.push_back(exp_add(weight+ratio_new,normlold+ratio_old)-rationorm);
		}
	}

	// partio[k]=log[1/n_k]
	std::vector<double> pratio(nreplica,0);
	// 1/n_k=n_1*\prod_j^k(m_j)
	for(unsigned i=0;i!=nreplica-1;++i)
		pratio[i+1]=pratio[i]+ratio[i];

	// As in the paper m_k=n_{k+1}/n_k
	// So here weights[k]=log[n_k]=log{1/pratio[k]}=-log{pratio[k]}
	
	//~ double rescalefb=logN+fb[rctid]-betak[rctid]*peshift;

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
				odebug.printField("ratio",int_ratios[i]);
			else
				odebug.printField("temperature",int_temps[i]);
			odebug.printField("betak",betak[i]);
			odebug.printField("fb",fb[i]);
			odebug.printField("old_fb",old_fb[i]);
			odebug.printField("rbfb",rbfb[i]);
			odebug.printField("rct",fb_rct[i]);
			odebug.printField("rct_rbfb",(fb_rct[i]-old_fb[i])/betak[i]);
			odebug.printField();
		}
		odebug.printf("--- FB update ---\n");
	}
	if(is_nan)
		plumed_merror("FB value become NaN. See debug file or \"error.data\" file for detailed information.");
}

void ITS_Bias::fb_partition()
{
	for(unsigned i=0;i!=nreplica;++i)
		rbfb[i]-=std::log(static_cast<double>(update_step));
		
	//~ if(is_first)
	//~ {
		//~ double rbfbsum=rbfb[0];
		//~ zbe[0]=fb_ratios[0]-fb[0];
		//~ double zisum=zbe[0];
		//~ for(unsigned i=1;i!=nreplica;++i)
		//~ {
			//~ exp_added(rbfbsum,rbfb[i]);
			//~ zbe[i]=fb_ratios[i]-fb[i];
			//~ exp_added(zisum,zbe[i]);
		//~ }
		//~ for(unsigned i=0;i!=nreplica;++i)
			//~ zbe[i]+=rbfbsum-zisum;
	//~ }

	if(mcycle==start_cycle)
	{
		for(unsigned i=0;i!=nreplica;++i)
			zbe[i]=rbfb[i];
	}
	else
	{
		if(!use_ema||mcycle<ema_days)
		{
			double ww=1.0/mcycle;
			double wt=std::log(ww);
			double wy=std::log(1.0-ww);
			for(unsigned i=0;i!=nreplica;++i)
				zbe[i]=exp_add(zbe[i]+wy,rbfb[i]+wt);
		}
		else
		{
			for(unsigned i=0;i!=nreplica;++i)
				zbe[i]=exp_add(zbe[i]+wemay,rbfb[i]+wemat);
		}
	}

	double zbesum=zbe[0];
	for(unsigned i=1;i!=nreplica;++i)
		exp_added(zbesum,zbe[i]);

	std::vector<double> pb(nreplica);
	for(unsigned i=0;i!=nreplica;++i)
		pb[i]=zbe[i]-zbesum;
	
	for(unsigned i=0;i!=nreplica;++i)
		fb[i]=fb_ratios[i]-pb[i];

	double fb0=fb[0];
	for(unsigned i=0;i!=nreplica;++i)
		fb[i]-=fb0;
}

void ITS_Bias::zbu_iteration()
{
	for(unsigned i=0;i!=nreplica;++i)
		rbzb[i]-=std::log(static_cast<double>(update_step));

	if(mcycle==start_cycle&&!read_norm)
	{
		for(unsigned i=0;i!=nreplica;++i)
			zbu[i]=rbzb[i];
	}
	else
	{
		if(mcycle<ema_days)
		{
			double ww=1.0/mcycle;
			double wt=std::log(ww);
			double wy=std::log(1.0-ww);
			for(unsigned i=0;i!=nreplica;++i)
				zbu[i]=exp_add(zbu[i]+wy,rbzb[i]+wt);
		}
		else
		{
			for(unsigned i=0;i!=nreplica;++i)
				zbu[i]=exp_add(zbu[i]+wemay,rbzb[i]+wemat);
		}
	}
}

void ITS_Bias::output_fb()
{
	if(mcycle%fb_stride==0)
	{
		if(abs_ratio)
			ofb.printField("ITERATE_METHOD","ABSOLUTE_RATIO");
		else
			ofb.printField("ITERATE_METHOD","RELATIVE_RATIO");
		ofb.printField("BOLTZMANN_CONSTANT",kB);
		ofb.printField("ITERATE_STEP",int(mcycle));
		ofb.printField("PACE",int(update_step));
		ofb.printField("PESHIFT",peshift);
		ofb.printField("COEFFICIENT_TYPE",(is_set_temps?"TEMPERATURE":"RATIO"));
		ofb.printField("NREPLICA",int(nreplica));
		ofb.printField("TEMP_RATIO_ENERGY",temp_ratio_energy);

		for(unsigned i=0;i!=nreplica;++i)
		{
			std::string id;
			Tools::convert(i,id);
			ofb.printField("index",id);
			if(is_set_ratios)
				ofb.printField("ratio",int_ratios[i]);
			else
				ofb.printField("temperature",int_temps[i]);
			ofb.printField("fb_value",fb[i]);
			if(abs_ratio)
				ofb.printField("Z_b(E)",zbe[i]);
			else
				ofb.printField("norm_value",norml[i]);
			ofb.printField("Z_b(U)",zbu[i]);
			ofb.printField();
		}
		ofb.printf("#!-----END-OF-FB-COEFFICIENTS-----\n\n");
		ofb.flush();
	}
	if(fbtrj_output&&mcycle%fbtrj_stride==0)
	{
		ofbtrj.fmtField(" %f");
		ofbtrj.printField("step",int(mcycle));

		for(unsigned i=0;i!=nreplica;++i)
		{
			std::string id;
			Tools::convert(i,id);
			std::string fbid="FB"+id;
			
			if(fbtrj_output)
				ofbtrj.printField(fbid,fb[i]);
		}
		ofbtrj.printField();
		ofbtrj.flush();
	}
	if(rct_output&&mcycle%rct_stride==0)
	{
		orct.fmtField(" %f");
		orct.printField("step",int(mcycle));
		for(unsigned i=0;i!=nreplica;++i)
		{
			std::string id;
			Tools::convert(i,id);
			std::string fbid="RCT"+id;
			if(rct_output)
				orct.printField(fbid,fb_rct[i]);
		}
		orct.printField();
		orct.flush();
	}

	if(bias_output&&mcycle%bias_stride==0)
		output_bias();
}

double ITS_Bias::find_rw_id(double rwtemp)
{
	if(rwtemp<templ)
		plumed_merror("the reweight temperature is lower than the minimal temperature of ITS");
	if(rwtemp>temph)
		plumed_merror("the reweight temperature is larger than the maximal temperature of ITS");
	
	double dtl(0),dth(1);
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
	
	return double(rwid)+dtl/(dtl+dth);
}

double ITS_Bias::find_rw_fb(double real_id)
{	
	unsigned rwid=unsigned(real_id);
	
	double fbl=fb[rwid];
	double fbh=fb[rwid+1];
	double dfb=fbh-fbl;
	
	double rtl=real_id-rwid;
	double rth=rwid+1.0-real_id;
	
	double rwfb;
	if(rtl<rth)
		rwfb=fbl+dfb*rtl;
	else
		rwfb=fbh-dfb*rth;
	return rwfb;
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
    ofile.setHeavyFlush();
}

void ITS_Bias::setupOFiles(std::vector<std::string>& fnames, std::vector<OFile*>& OFiles, const bool multi_sim_single_files)
{
	OFiles.resize(fnames.size(),NULL);
	for(unsigned i=0; i!=fnames.size();++i)
	{
		OFiles[i] = new OFile();
		OFiles[i]->link(*this);
		if(multi_sim_single_files)
		{
			unsigned int r=0;
			if(comm.Get_rank()==0)
				{r=multi_sim_comm.Get_rank();}
			comm.Bcast(r,0);
			if(r>0)
				{fnames[i]="/dev/null";}
			OFiles[i]->enforceSuffix("");
		}
		OFiles[i]->open(fnames[i]);
		OFiles[i]->setHeavyFlush();
	}
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
	
	unsigned _ntemp=nreplica;

	unsigned read_count=0;
	
	std::string iter_method="RELATIVE_RATIO";
	
	while(ifb)
	{
		if(ifb.FieldExist("ITERATE_METHOD"))
			ifb.scanField("ITERATE_METHOD",iter_method);

		if(ifb.FieldExist("BOLTZMANN_CONSTANT"))
			ifb.scanField("BOLTZMANN_CONSTANT",_kB);

		if(ifb.FieldExist("ITERATE_STEP"))
		{
			int tmc=0;
			ifb.scanField("ITERATE_STEP",tmc);
			mcycle=tmc;
		}
		
		if(iter_method=="RELATIVE_RATIO"||iter_method=="Traditional")
		{
			if(abs_ratio)
				mcycle=0;
		}
		else if(iter_method=="ABSOLUTE_RATIO"||iter_method=="Partition")
		{
			if(!abs_ratio)
				mcycle=0;
		}
		else
			plumed_merror("unknown iterate method: \""+iter_method+"\"");
		
		if(ifb.FieldExist("PACE"))
		{
			int tmc=0;
			ifb.scanField("PACE",tmc);
			update_step=tmc;
		}

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
		
		if(ifb.FieldExist("TEMP_RATIO_ENERGY"))
			ifb.scanField("TEMP_RATIO_ENERGY",temp_ratio_energy);
		
		int_ratios.resize(_ntemp);
		int_temps.resize(_ntemp);
		fb.resize(_ntemp);
		
		if(ifb)
		{
			if(!ifb.FieldExist("fb_value"))
				plumed_merror("cannot found \"fb_value\" in file \""+fb_input+"\"");
			if(is_set_temps&&!ifb.FieldExist("temperature"))
				plumed_merror("cannot found \"temperature\" in file \""+fb_input+"\"");
			if(is_set_ratios&&!ifb.FieldExist("ratio"))
				plumed_merror("cannot found \"ratio\" in file \""+fb_input+"\"");
		}
		else
			break;

		if(!is_const)
		{
			norml.assign(_ntemp,1e-38);
			zbe.assign(_ntemp,1e-38);
			if(mcycle>0&&(ifb.FieldExist("norm_value")||ifb.FieldExist("Z_b(E)")))
				read_norm=true;
		}
		
		if(ifb.FieldExist("Z_b(U)"))
		{
			zbu.resize(_ntemp);
			read_zbu=true;
		}
		
		double tmpfb;
		
		for(unsigned i=0;i!=_ntemp;++i)
		{
			ifb.scanField("fb_value",tmpfb);
			double tmpt;
			if(is_set_ratios)
			{
				ifb.scanField("ratio",tmpt);
				int_ratios[i]=tmpt;
				int_temps[i]=sim_temp*tmpt;
			}
			else
			{
				ifb.scanField("temperature",tmpt);
				int_temps[i]=tmpt;
				int_ratios[i]=tmpt/sim_temp;
			}
			
			if(read_norm)
			{
				double tmpnb;
				if(abs_ratio)
				{
					ifb.scanField("Z_b(E)",tmpnb);
					zbe[i]=tmpnb;
				}
				else
				{
					ifb.scanField("norm_value",tmpnb);
					norml[i]=tmpnb;
				}
			}
			
			if(read_zbu)
			{
				double zk;
				ifb.scanField("Z_b(U)",zk);
				zbu[i]=zk;
			}
			fb[i]=tmpfb;
			
			ifb.scanField();
		}
		++read_count;
	}
	
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
		temph=sim_temp*ratiol;
		templ=sim_temp*ratioh;
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
		ratioh=templ/sim_temp;
		ratiol=temph/sim_temp;
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
	peshift=new_shift;
	set_peshift_ratio();
}

void ITS_Bias::set_peshift_ratio()
{
	for(unsigned i=0;i!=nreplica-1;++i)
		peshift_ratio[i]=(betak[i]-betak[i+1])*peshift;
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

}
}
