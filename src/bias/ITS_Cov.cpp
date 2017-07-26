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
#include "Bias.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <ctime>
#include <cmath>

using namespace std;


namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS ITS_VES 
/*

Code by Y. Isaac Yang, Ph.D.
ETH Zurich, USI-Campus
Email: yi.yang@phys.chem.ethz.ch

Last modified: October 19, 2016
*/
//+ENDPLUMEDOC

class ITS_VES: public Bias
{
private:
	const unsigned update_start;
	unsigned ntemp;
	double templ;
	double temph;
	double kB;
	double kBT;
	double beta0;
	double sim_temp;
	double peshift;
	double shift_energy;

	std::vector<double> int_temps;
	std::vector<double> betak;

	std::vector<double> fb;
	std::vector<double> norml;
	
	std::vector<double> gU;
	std::vector<double> gf;
	std::vector<double> bgf;
	std::vector<double> rbfb;

	std::vector<double> energy_disperse;
	std::vector<double> energy_prob;
	std::vector<double> fb_iter;

	std::vector<double> rw_temp;
	std::vector<double> rw_beta;
	std::vector<double> rw_factor;
	std::vector<std::string> rw_tstr;

	unsigned step;
	unsigned norm_step;
	unsigned mcycle;
	unsigned tar_id;
	double fb_init;
	double fb_bias;
	double rb_fac1;
	double rb_fac2;
	double rb_slow;
	double rb_slow2;
	double bias_min;
	double bias_max;

	bool is_const;
	bool is_output;
	bool read_fb;
	bool read_norm;
	bool read_iter;
	bool is_ves;
	bool only_1st;
	bool use_mw;
	bool bias_output;
	bool is_debug;
	unsigned update_step;
	unsigned output_step;
	unsigned output_start;
	unsigned start_cycle;
	unsigned fb_stride;
	unsigned fbtrj_stride;
	unsigned bias_stride;
	unsigned bias_bins;
	unsigned debug_step;

	std::string fb_file;
	std::string fb_trj;
	std::string norm_trj;
	std::string peshift_trj;
	std::string rw_file;
	std::string fb_input;
	std::string bias_file;
	std::string debug_file;

	OFile ofb;
	OFile ofbtrj;
	OFile onormtrj;
	OFile opstrj;
	OFile ofw;
	OFile obias;
	OFile odebug;

	double gfsum;
	double bgfsum;
	double min_ener;

	unsigned ntarget;
	double tar_min;
	double tar_max;
	double tar_ratio;
	double var_step_size;
	double dEt;
	double d_pot;

	Value* valueEff;
	Value* valueForce;
	Value* valueReweight;

	inline void update_rbfb();
	void iterate_fb();

	void fb_variational();
	void output_fb();
	void mw_merge_rbfb();
	void output_bias();

	inline double force_scale(double v2);
	inline void coe_rescale(double shift,std::vector<double>& coe);
	inline void fb_rescale(double shift){coe_rescale(shift,fb);}
	inline void iter_rescale(double shift){coe_rescale(shift,fb_iter);}
	inline void change_peshift(double new_shift);

	inline double exp_add(double exp1,double exp2);
	inline void exp_added(double& expsum,double expvlaue);

	inline double exp_minus(double exp1,double exp2,bool& sign);
	inline double exp_calc(double exp1,double exp2,bool sign1,bool sign2,bool& sign_fin);

	inline std::vector<double> omega2_alpha(const std::vector<std::vector<double> >& hessian,const std::vector<double>& dalpha);
	inline std::vector<double> omega2_alpha(const std::vector<double>& hessian,const std::vector<double>& dalpha);

public:
	explicit ITS_VES(const ActionOptions&);
	~ITS_VES();
	void calculate();
	inline double reweight(double rw_beta,double energy);
	double reweight(double rw_beta){return exp(-rw_beta*shift_energy-gfsum);}
	static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ITS_VES,"ITS_VES")

void ITS_VES::registerKeywords(Keywords& keys)
{
	Bias::registerKeywords(keys);
	keys.addOutputComponent("effective","default","the instantaneous value of the effective potential");
	keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
	keys.addOutputComponent("force","default","the instantaneous value of the bias force");
	keys.addOutputComponent("reweight","default","instantaneous value of the default reweighting factor");
	ActionWithValue::useCustomisableComponents(keys);
	keys.use("ARG");

	keys.add("compulsory","NTEMP","100","the number of the temperatures");
	keys.add("compulsory","PESHIFT","0.0","the shift value of potential energy");
	keys.add("compulsory","TEMP_MIN","273.0","the lower bounds of the temperatures");
	keys.add("compulsory","TEMP_MAX","373.0","the upper bounds of the temperatures");
	keys.add("compulsory","PACE","1000","the frequency for updating fb value");

	keys.add("compulsory","START_CYCLE","0","the start step for fb updating");
	keys.add("compulsory","FB_INIT","-0.01","the default value for fb initializing");
	keys.add("compulsory","FB_BIAS","0.0","the bias value for fb updating");
	keys.add("compulsory","RB_FAC1","0.5","the ratio of the average value of rb");
	keys.add("compulsory","RB_FAC2","0.0","the ratio of the old steps in rb updating");
	keys.add("compulsory","RB_SLOW","1.0","the constant that make the iteration of fb slow down");

	keys.add("compulsory","FB_TRAJ","fbtrj.data","a file recording the evolution of fb values");
	keys.add("compulsory","NORM_TRAJ","normltrj.data","a file recording the evolution of normalize factors");
	keys.add("compulsory","PESHIFT_TRAJ","peshift_trj.data","a file recording the evolution of peshift");
	keys.add("compulsory","FB_OUT","fb.data","a file to record the new fb values when they are update");
	keys.add("compulsory","FB_OUT_STRIDE","1","the frequency to output the new fb values");
	keys.add("compulsory","FBTRAJ_STRIDE","1","the frequency to record the evoluation of fb values");

	keys.addFlag("FB_FIXED",false,"to fix the fb value without update");
	keys.addFlag("MULTIPLE_WALKERS",false,"use multiple walkers");
	keys.addFlag("VARIATIONAL",false,"use variational approach to iterate the fb value");
	keys.addFlag("ONLY_FIRST_ORDER",false,"only to calculate the first order of Omega (gradient) in variational approach");

	keys.add("optional","SIM_TEMP","the temperature used in the simulation");
	keys.add("optional","RW_TEMP","the temperatures used in the calcaulation of reweighting factors");
	keys.add("optional","RW_FILE","a file of the reweighting factors");
	keys.add("optional","RW_STRIDE","the frequency of reweighting factor output");
	keys.add("optional","OUTPUT_START","the start step of reweighting factor output");
	keys.add("optional","FB_READ_FILE","a file of reading fb values (include temperatures and peshift)");
	keys.add("optional","BIAS_FILE","a file to output the function of bias potential");
	keys.add("optional","BIAS_STRIDE","the frequency to output the bias potential");
	keys.add("optional","BIAS_MIN","the minimal value of coordinate of the bias potential function");
	keys.add("optional","BIAS_MAX","the maximal value of coordinate of the bias potential function");
	keys.add("optional","BIAS_BIN","the number of bins to record the bias potential function");
	keys.add("optional","DEBUG_FILE","a file of output debug information");
	keys.add("optional","DEBUG_STEP","the frequency for debug file output");

	keys.add("optional","TARGET_MIN","the minimal potential energy of the uniform target distribution used in the variational approach");
	keys.add("optional","TARGET_MAX","the maximal potential energy of the uniform target distribution used in the variational approach");
	keys.add("optional","PROB_RATIO","the ratio between the probability of maximal and minimal energys");
	keys.add("optional","TARGET_NUMBER","the number of segment of the uniform target distribution used in the variational approach");
	keys.add("optional","VAR_STEP_SIZE","the step size of iteration used in the variational approach");
}

ITS_VES::~ITS_VES()
{
	if(!is_const)
	{
		ofb.close();
		ofbtrj.close();
		onormtrj.close();
		opstrj.close();
	}
	if(is_output)
		ofw.close();
	if(is_debug)
		odebug.close();
}

ITS_VES::ITS_VES(const ActionOptions& ao):
	PLUMED_BIAS_INIT(ao),update_start(0),step(0),norm_step(0),mcycle(0),
	is_const(false),is_output(false),
	read_fb(false),read_norm(false),read_iter(false),
	is_ves(false),only_1st(false),bias_output(false),is_debug(false),
	output_start(0),bias_stride(1),debug_step(1),min_ener(1e38),tar_ratio(1.0),var_step_size(1.0)
{
	if(getNumberOfArguments()>1)
		error("this edition of ITS can only accept one CV");

	kB=plumed.getAtoms().getKBoltzmann();

	parse("PESHIFT",peshift);
	parseFlag("VARIATIONAL",is_ves);
	parseFlag("FB_FIXED",is_const);

	parse("FB_READ_FILE",fb_input);
	
	double _ntemp,_templ,_temph;
	parse("NTEMP",_ntemp);
	parse("TEMP_MIN",_templ);
	parse("TEMP_MAX",_temph);
	parse("START_CYCLE",start_cycle);
	
	double _kB=kB;
	double _peshift=peshift;
	if(fb_input.size()>0)
	{
		read_fb=true;
		IFile ifb;
		ifb.open(fb_input);
		ifb.allowIgnoredFields();
		
		string iter_method="Tradtional";
		if(ifb.FieldExist("ITERATE_METHOD"))
			ifb.scanField("ITERATE_METHOD",iter_method);
			
		if(iter_method=="Traditional")
		{
			if(is_ves)
				error("The iterate method is variational but the fb input file is a traditional input file");
		}
		else if(iter_method=="Variational")
		{
			if(!is_ves)
				error("The iterate method is traditional but the fb input file is a variational input file");
		}
		else
			error("unknown iterate method: \""+iter_method+"\"");

		
		if(ifb.FieldExist("ITERATE_STEP"))
		{
			int tmc;
			ifb.scanField("ITERATE_STEP",tmc);
			mcycle=tmc;
		}
		
		if(ifb.FieldExist("BOLTZMANN_CONSTANT"))
			ifb.scanField("BOLTZMANN_CONSTANT",_kB);
		
		if(ifb.FieldExist("PESHIFT"))
			ifb.scanField("PESHIFT",_peshift);
		peshift=_peshift;
		if(fabs(_kB/kB-1)>1.0e-8)
			peshift*=kB/_kB;

		if(!ifb.FieldExist("fb_value"))
			error("cannot found \"fb_value\" in file \""+fb_input+"\"");
		if(!ifb.FieldExist("temperature"))
			error("cannot found \"temperature\" in file \""+fb_input+"\"");
		double tmpfb;
		int_temps.resize(0);

		if(!is_ves&&!is_const&&ifb.FieldExist("norm_value"))
			read_norm=true;
		if(is_ves&&ifb.FieldExist("iterate_fb"))
			read_iter=true;
		
		while(ifb.scanField("fb_value",tmpfb))
		{
			double tmpt;
			ifb.scanField("temperature",tmpt);
			int_temps.push_back(tmpt);
			
			if(!is_ves&&read_norm)
			{
				double tmpnb;
				ifb.scanField("norm_value",tmpnb);
				norml.push_back(tmpnb);
			}
			if(is_ves&&read_iter)
			{
				double tmpif;
				ifb.scanField("iterate_fb",tmpif);
				fb_iter.push_back(tmpif);
			}
			
			fb.push_back(tmpfb);
			ifb.scanField();
		}
		
		ntemp=int_temps.size();
		templ=temph=int_temps[0];
		for(unsigned i=1;i!=int_temps.size();++i)
		{
			if(int_temps[i]>temph)
				temph=int_temps[i];
			if(int_temps[i]<templ)
				templ=int_temps[i];
		}
		ifb.close();
	}
	else
	{
		if(_ntemp<1)
			error("the number of temperatures much be larger than 1");
		ntemp=_ntemp;
		
		if(_templ<0)
			error("TEMP_MIN must be large than 0");
		templ=_templ;
		
		if(_temph<0)
			error("TEMP_MAX must be large than 0");
		temph=_temph;
		if(templ>=temph)
			error("the value of TEMP_MAX must be large than TEMP_MIN");
		double dt=(temph-templ)/(ntemp-1);
		for(unsigned i=0;i!=ntemp;++i)
			int_temps.push_back(templ+i*dt);
		mcycle=start_cycle;
	}

	parse("SIM_TEMP",sim_temp);
	if(sim_temp>0)
	{
		kBT=kB*sim_temp;
	}
	else
	{
		kBT=plumed.getAtoms().getKbT();
		sim_temp=kBT/kB;
	}
	beta0=1.0/kBT;

	if(sim_temp>=temph||sim_temp<=templ)
	{
		string sts,stl,sth;
		Tools::convert(sim_temp,sts);
		Tools::convert(templ,stl);
		Tools::convert(temph,sth);
		error("the value of SIM_TEMP("+sts+") must between TEMP_MIN("+stl+") and TEMP_MAX("+sth+")");
	}

	for(unsigned i=0;i!=ntemp;++i)
	{
		betak.push_back(1.0/(kB*int_temps[i]));
		if(!read_fb)
			fb.push_back(-exp(fb_init*ntemp));
	}
	if(!read_norm)
		norml.assign(ntemp,0);

	gU.assign(ntemp,0);
	gf.assign(ntemp,0);
	bgf.assign(ntemp,0);
	rbfb.assign(ntemp,0);

	parse("PACE",update_step);
	if(update_step==0)
		error("PACE cannot be 0");

	parse("FB_INIT",fb_init);
	parse("FB_BIAS",fb_bias);
	parse("RB_FAC1",rb_fac1);
	parse("RB_FAC2",rb_fac2);
	parse("RB_SLOW",rb_slow);
	if(rb_slow<=0)
		error("RB_SLOW must be larger than 0!");
	rb_slow2=std::log(rb_slow);
	parse("START_CYCLE",start_cycle);

	parse("FB_TRAJ",fb_trj);
	parse("NORM_TRAJ",norm_trj);
	parse("PESHIFT_TRAJ",peshift_trj);
	parse("FB_OUT",fb_file);
	parse("FB_OUT_STRIDE",fb_stride);
	parse("FBTRAJ_STRIDE",fbtrj_stride);

	if(!is_const)
	{
		ofb.link(*this);
		ofb.open(fb_file);

		ofbtrj.link(*this);
		ofbtrj.open(fb_trj);

		onormtrj.link(*this);
		onormtrj.open(norm_trj);
		
		opstrj.link(*this);
		opstrj.open(peshift_trj);
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
			error("BIAS_MAX must be larger than bias_min");
		parse("BIAS_BIN",bias_bins);
		if(bias_bins<=1)
			error("BIAS_BIN must be larger than 1");
		parse("BIAS_STRIDE",bias_stride);
		d_pot=(bias_max-bias_min)/(bias_bins-1);
	}
	
	parseFlag("MULTIPLE_WALKERS",use_mw);

	parse("TARGET_MIN",tar_min);
	parse("TARGET_MAX",tar_max);
	parse("PROB_RATIO",tar_ratio);
	parse("TARGET_NUMBER",ntarget);
	parse("VAR_STEP_SIZE",var_step_size);
	parseFlag("ONLY_FIRST_ORDER",only_1st);

	if(is_ves)
	{	
		if(!read_iter)
		{
			double fb_zero=fb[0];
			for(unsigned i=0;i!=ntemp;++i)
			{
				fb[i]-=fb_zero;
				fb_iter.push_back(fb[i]);
			}
		}

		dEt=(tar_max-tar_min)/(ntarget-1);

		double slope=(tar_ratio-1.0)/(ntarget-1);
		double probsum=(1.0+tar_ratio)/2*ntarget;
		for(unsigned i=0;i!=ntarget;++i)
		{
			energy_disperse.push_back(tar_min+i*dEt);
			energy_prob.push_back((1.0+i*slope)/probsum);
		}
	}

	parseVector("RW_TEMP",rw_temp);
	parse("RW_FILE",rw_file);
	parse("RW_STRIDE",output_step);
	parse("OUTPUT_START",output_start);

	if(rw_file.size()==0&&rw_temp.size()>0)
		error("RW_TEMP is setted but RW_FILE is not setted");

	if(rw_file.size()>0)
	{
		is_output=true;
		if(rw_temp.size()>0)
		{
			for(unsigned i=0;i!=rw_temp.size();++i)
			{
				if(rw_temp[i]<=temph&&rw_temp[i]>=templ)
					rw_beta.push_back(1.0/(kB*rw_temp[i]));
				else
					error("the reweighting temperature must between TEMP_MIN and TEMP_MAX");
			}
		}
		else
		{
			rw_temp.push_back(sim_temp);
			rw_beta.push_back(beta0);
		}
		rw_factor.resize(rw_temp.size());

		ofw.link(*this);
		ofw.open(rw_file);
		if(!is_const)
			ofw.printf("#! Warning: As the fb values have not been fixed yet, the reweighting factor may be meaningless!\n");
		for(unsigned i=0;i!=rw_temp.size();++i)
		{
			string tt;
			Tools::convert(i,tt);
			tt="rw_temperature_"+tt;
			ofw.addConstantField(tt).printField(tt,rw_temp[i]);
		}

		if(output_step==0)
			output_step=update_step;
	}

	parse("DEBUG_FILE",debug_file);
	parse("DEBUG_STEP",debug_step);
	if(debug_file.size()>0)
	{
		is_debug=true;
		odebug.link(*this);
		odebug.open(debug_file);
	}

	checkRead();

	addComponent("effective"); componentIsNotPeriodic("effective");
	valueEff=getPntrToComponent("effective");
	addComponent("force"); componentIsNotPeriodic("force");
	valueForce=getPntrToComponent("force");
	addComponent("reweight"); componentIsNotPeriodic("reweight");
	valueReweight=getPntrToComponent("reweight");

	log.printf("  Simulation Temperature: %f\n",sim_temp);
	log.printf("  Boltzmann Constant: %f\n",kB);
	log.printf("  Beta (1/kT): %f\n",beta0);

	if(read_fb)
		log.printf("  FB value read from file: %s\n",fb_input.c_str());
	log.printf("  Temperatues from %f to %f\n",templ,temph);
	log.printf("  Temperatue and fb value:\n");
	for(unsigned i=0;i!=ntemp;++i)
		log.printf("    %d\t%f\t%f\n",i,int_temps[i],fb[i]);
	log.printf("  Number of temperature segment: %d\n",ntemp);
	log.printf("  PESHIFT value: %f\n",peshift);
	if(fabs(_kB/kB-1)>1.0e-8)
		log.printf("    Original PESHIFT value: %f (with boltzmann constant %f)\n",_peshift,_kB);
	if(read_norm)
	{
		log.printf("  Normalization value read from file: %s\n",fb_input.c_str());
		for(unsigned i=0;i!=norml.size();++i)
			log.printf("    %d\t%f\n",i,norml[i]);
	}
	if(!is_const)
	{
		log.printf("  Using ITS_VES update\n");
		log.printf("    Frequence of FB value update: %d\n",update_step);
		log.printf("    FB output file: %s\n",fb_file.c_str());
		log.printf("    FB trajectory file: %s\n",fb_trj.c_str());
		log.printf("    Normalized factors trajectory file: %s\n",norm_trj.c_str());
		log.printf("    Potential energy shift trajectory file: %s\n",peshift_trj.c_str());
		log.printf("\n");
	}
	if(bias_output)
	{
		log.printf("  Output bias potential function\n");
		log.printf("    Bias output file: %s\n",bias_file.c_str());
		log.printf("    Output range from %f to %f\n",bias_min,bias_max);
		log.printf("    Output bins: %d\n",bias_bins);
		log.printf("    Frequence of output: %d\n",bias_stride);
		log.printf("\n");
	}
	if(use_mw)
	{
		log.printf("  Using multiple walkers");
		log.printf("   number of walkers: %d\n",multi_sim_comm.Get_size());
		log.printf("   walker number: %d\n",multi_sim_comm.Get_rank());
		log.printf("\n");
	}
	if(is_output)
	{
		log.printf("  Reweighting factor output file: %s\n",rw_file.c_str());
		log.printf("  Reweighting factor at temperature:");
		for(unsigned i=0;i!=rw_temp.size();++i)
			log.printf(" %f",rw_temp[i]);
		log.printf("\n");
		log.printf("  Frequence of reweighting factor output: %d\n",output_step);
	}
	if(is_debug)
		log.printf("  Debug Mod is open, output file: %s\n",debug_file.c_str());

	log<<"Bibliography "<<
		plumed.cite("Y. Q. Gao, J. Chem. Phys. 128, 064105 (2008)")<<
		plumed.cite("Y. Q. Gao, J. Chem. Phys. 128, 134111 (2008)");
	log<<"\n";
}

void ITS_VES::calculate()
{
	if(!is_const && step%update_step==update_start && step>update_start)
	{
		++mcycle;
		if(is_ves)
			fb_variational();
		else
			iterate_fb();
		norm_step=0;
	}
	if(bias_output&&step==0)
		output_bias();

	double energy=getArgument(0);
	shift_energy=energy+peshift;

	if(energy<min_ener)
		min_ener=energy;

	for(unsigned i=0;i!=ntemp;++i)
	{
		gU[i]=-betak[i]*shift_energy;
		gf[i]=gU[i]+fb[i];
		bgf[i]=gf[i]+std::log(betak[i]);
	}
	gfsum=gf[0];
	bgfsum=bgf[0];
	for(unsigned i=1;i!=ntemp;++i)
	{
		exp_added(gfsum,gf[i]);
		exp_added(bgfsum,bgf[i]);
	}

	double eff_pot=-gfsum/beta0;
	double bias=eff_pot-shift_energy;
	//~ double force=exp(bgfsum-gfsum)/beta0-1;
	double force=1-exp(bgfsum-gfsum)/beta0;
	double def_rwf=exp(-beta0*shift_energy-gfsum);

	if(is_output && step%output_step==output_start)
	{
		ofw.printField("time",getTimeStep()*step);
		ofw.printField("energy",energy);
		ofw.printField("peshift",peshift);
		for(unsigned i=0;i!=rw_factor.size();++i)
		{
			string tt;
			Tools::convert(i,tt);
			tt="rw_factor_"+tt;
			rw_factor[i]=reweight(rw_beta[i]);
			ofw.printField(tt,rw_factor[i]);
		}
		ofw.printField();
		ofw.flush();
	}

	if(is_debug&&step%debug_step==update_start)
	{
		odebug.addConstantField("STEP").printField("STEP",int(step));
		odebug.addConstantField("ENERGY").printField("ENERGY",energy);
		odebug.addConstantField("PESHIFT").printField("PESHIFT",peshift);
		odebug.addConstantField("gfsum").printField("gfsum",gfsum);
		odebug.addConstantField("bgfsum").printField("bgfsum",bgfsum);
	}

	if(!is_const)
	{
		if(!is_ves||(energy>=tar_min&&energy<=tar_max))
		{
			update_rbfb();
			++norm_step;
		}
	}

	if(is_debug)
		odebug.flush();

	setBias(bias);
	setOutputForce(0,force);

	valueEff->set(eff_pot);
	valueForce->set(force);
	valueReweight->set(def_rwf);

	++step;
}

inline void ITS_VES::update_rbfb()
{
	const unsigned& stride=comm.Get_size();
	const unsigned& rank=comm.Get_rank();
	// In traditional iterate process:
	// rbfb[k]=log[\sum_t(p_k)]; p_k=P_k/[\sum_k(P_k)]; P_k=int_
	// the summation record the data of each the update steps (default=100)
	// In variational iterate process:
	// rbfb[i]=log[(\sum_t(\beta*(\partial V_bias(U;a)/\partial a_i)))_V(a)]
	//~ if(step%update_step==update_start)
	if(norm_step==0)
	{
		for(unsigned i=rank;i<ntemp;i+=stride)
			rbfb[i]=gf[i]-gfsum;
	}
	else
	{
		for(unsigned i=rank;i<ntemp;i+=stride)
			exp_added(rbfb[i],gf[i]-gfsum);
	}
	if(is_debug&&step%debug_step==update_start)
	{
		for(unsigned i=0;i!=ntemp;++i)
		{
			string id;
			Tools::convert(i,id);
			odebug.printField("index",id);
			odebug.printField("stride",int(stride));
			odebug.printField("temperature",int_temps[i]);
			odebug.printField("fb",fb[i]);
			odebug.printField("gf",gf[i]);
			odebug.printField("bgf",bgf[i]);
			odebug.printField("rbfb",rbfb[i]);
			odebug.printField();
		}
	}
}

// Y. Q. Gao, J. Chem. Phys. 128, 134111 (2008)
void ITS_VES::iterate_fb()
{
	comm.Sum(rbfb);
	if(use_mw)
		mw_merge_rbfb();
	if(-1.0*min_ener>peshift)
		change_peshift(-1.0*min_ener);
	for(unsigned i=0;i!=ntemp;++i)
		rbfb[i]-=rb_slow2;

	std::vector<double> rb;
	std::vector<double> ratio1;

	for(unsigned i=0;i!=ntemp-1;++i)
	{
		// rb[k]=log[f(p_k,p_{k+1})]
		rb.push_back((rbfb[i]+rbfb[i+1])*rb_fac1+mcycle*rb_fac2);
		// ratio_old[k]=log[m_k(t-1)], m_k=n_k/n_{k+1}
		// (Notice that in the paper m_k=n_{k+1}/n_k)
		ratio1.push_back(fb[i]-fb[i+1]);
	}

	// normal=log[W_k(t)], normalold=log[W_k(t-1)]
	std::vector<double> normlold(norml);
	if(mcycle==1&&!read_norm)
	{
		for(unsigned i=0;i!=ntemp-1;++i)
		{
			norml[i]=rb[i];
			normlold[i]=-1e38;
		}
	}
	else
	{
		// W_k(t)=\sum_t[f(p_k,p_{k+1})]
		// the summation record the data of ALL the update steps
		for(unsigned i=0;i!=ntemp-1;++i)
			exp_added(norml[i],rb[i]);
	}
	
	// ratio[k]=log[m_k(t)]
	std::vector<double> ratio;
	// m_k(t)=[c_bias*f(p_k,p_{k+1})*p_k/p_{k+1}+m_k(t-1)*W_k(t-1)]/W_k(t)
	for(unsigned i=0;i!=ntemp-1;++i)
		ratio.push_back(ratio1[i]-norml[i]+exp_add(normlold[i],rbfb[i+1]-rbfb[i]+rb[i]+fb_bias));

	// partio[k]=log[1/n_k]
	std::vector<double> pratio(ntemp,0);
	// 1/n_k=n_1*\prod_j^k(m_j)
	for(unsigned i=0;i!=ntemp-1;++i)
		pratio[i+1]=pratio[i]+ratio[i];

	// As in the paper m_k=n_{k+1}/n_k
	// So here weights[k]=log[n_k]=log{1/pratio[k]}=-log{pratio[k]}
	for(unsigned i=0;i!=ntemp;++i)
	{
		fb[i]=-pratio[i];
		rbfb[i]=0;
	}
	output_fb();
}

void ITS_VES::fb_variational()
{
	comm.Sum(rbfb);
	if(use_mw)
		mw_merge_rbfb();
	if(-1.0*min_ener>peshift)
		change_peshift(-1.0*min_ener);

	// The average of 1st order drivative of bias potential in target distribution
	std::vector<double> deriv_ps(ntemp,0);
	// The average of 2nd order drivative of bias potential in target distribution
	std::vector<double> deriv2_ps(ntemp*(ntemp+1)/2,0);
	std::vector<std::vector<double> > deriv2_ps2(ntemp,std::vector<double>(ntemp,0));

	const unsigned& stride=comm.Get_size();
	const unsigned& rank=comm.Get_rank();
	for(unsigned idu=rank;idu<ntarget;idu+=stride)
	{
		double ener=energy_disperse[idu]+peshift;
		double prob=energy_prob[idu];

		double gfsum_;
		std::vector<double> gf_;

		gf_.push_back(-betak[0]*ener+fb[0]);
		gfsum_=gf_[0];
		for(unsigned i=1;i!=ntemp;++i)
		{
			double tmp=-betak[i]*ener+fb[i];
			gf_.push_back(tmp);
			exp_added(gfsum_,tmp);
		}

		std::vector<double>::iterator idv2=deriv2_ps.begin();
		for(unsigned i=0;i!=ntemp;++i)
		{
			deriv_ps[i]+=prob*exp(gf_[i]-gfsum_);
			if(!only_1st)
			{
				for(unsigned j=i;j!=ntemp;++j)
				{
					double dv2=prob*exp(gf_[i]+gf_[j]-2*gfsum_);
					(*(idv2++))+=dv2;
					deriv2_ps2[i][j]+=dv2;
					if(i!=j)
						deriv2_ps2[j][i]+=dv2;
				}
			}
		}
	}
	if(stride>1)
	{
		comm.Sum(deriv_ps);
		comm.Sum(deriv2_ps);
	}
	
	if(is_debug)
	{
		odebug.printf("--- driv2_ps2 ---\n");
		for(unsigned i=0;i!=ntemp;++i)
		{
			odebug.printField("temp_id",int(i));
			for(unsigned j=0;j!=ntemp;++j)
			{
				std::string title="deriv2_ps";
				string id;
				Tools::convert(j,id);
				title+=id;
				odebug.printField(title,deriv2_ps2[i][j]);
			}
			odebug.printField();
		}
		odebug.printf("--- END ---\n");
	}

	// The average of 1st order drivative of bias potential in simulation
	std::vector<double> deriv_aver(ntemp,0);
	// fb_diff[i]=fb_iter[i]-fb[i]=a(i)-<a(i)>
	std::vector<double> fb_diff(ntemp,0);
	// The gradient of Omega
	std::vector<double> gradient(ntemp,0);
	
	for(unsigned i=0;i!=ntemp;++i)
	{
		deriv_aver[i]=exp(rbfb[i])/norm_step;
		//~ deriv_ps[i]/=ntarget;
		gradient[i]=(deriv_aver[i]-deriv_ps[i])/beta0;
		if(!only_1st)
			fb_diff[i]=fb_iter[i]-fb[i];
	}

	// The hessian of Omega
	std::vector<std::vector<double> > hessian(ntemp,std::vector<double>(ntemp,0));
	if(!only_1st)
	{
		std::vector<double>::iterator idv2ps=deriv2_ps.begin();
		for(unsigned i=0;i!=ntemp;++i)
		{
			for(unsigned j=i;j!=ntemp;++j)
			{
				// <\partial V_bias(U;a)/\partial a_i>*<\partial V_bias(U;a)/\partial a_j>
				double deriv_aver_2=deriv_aver[i]*deriv_aver[j];
				hessian[i][j]=(*(idv2ps++)-deriv_aver_2)/beta0;
				if(i!=j)
					hessian[j][i]=hessian[i][j];
			}
		}
	}

	double fb_zero;
	vector<double> o2as;
	for(unsigned i=0;i!=ntemp;++i)
	{
		double fin_diff;
		if(only_1st)
			fin_diff=var_step_size*gradient[i];
		else
		{
			double o2a=0;
			for(unsigned j=0;j!=ntemp;++j)
				o2a+=hessian[i][j]*fb_diff[j];
			o2as.push_back(o2a);
			fin_diff=var_step_size*(gradient[i]+o2a);
		}

		if(i==0)
		{
			fb_zero=fin_diff;
			fin_diff=0;
		}
		else
			fin_diff-=fb_zero;

		fb_iter[i]-=fin_diff;
		fb[i]+=(fb_iter[i]-fb[i])/(mcycle+1);
		
		rbfb[i]=0;
	}
	
	if(is_debug)
	{
		odebug.printf("--- FB Variational ---\n");
		odebug.addConstantField("norm_step").printField("norm_step",int(norm_step));
		odebug.addConstantField("var_step_size").printField("var_step_size",var_step_size);
		odebug.addConstantField("peshift").printField("peshift",peshift);
		for(unsigned i=0;i!=ntemp;++i)
		{
			odebug.printField("index",int(i));
			odebug.printField("deriv_aver",deriv_aver[i]);
			odebug.printField("deriv_ps",deriv_ps[i]);
			odebug.printField("gradient",gradient[i]);
			odebug.printField("o2a",o2as[i]);
			odebug.printField("fb",fb[i]);
			odebug.printField("fb_iter",fb_iter[i]);
			odebug.printField();
		}
		odebug.printf("--- Variational END ---\n");
	}
	
	output_fb();
}

inline double ITS_VES::reweight(double rw_beta,double energy)
{
	energy+=peshift;
	double factor=0;
	for(unsigned i=0;i!=ntemp;++i)
		factor+=exp(fb[i]+(rw_beta-betak[i])*energy);
	factor=1.0/factor;
	return factor;
}

inline void ITS_VES::coe_rescale(double shift,std::vector<double>& coe)
{
	double coeshift=coe[0]+betak[0]*shift;
	for(unsigned i=0;i!=coe.size();++i)
		coe[i]+=betak[i]*shift-coeshift;
}

inline void ITS_VES::change_peshift(double new_shift)
{
	fb_rescale(new_shift-peshift);
	if(is_ves)
		iter_rescale(new_shift-peshift);
	peshift=new_shift;
}

inline double ITS_VES::exp_add(double exp1,double exp2)
{
	if(exp1>exp2)
		return exp1+std::log(1.0+exp(exp2-exp1));
	else
		return exp2+std::log(1.0+exp(exp1-exp2));
}

inline void ITS_VES::exp_added(double& expsum,double expvalue)
{
	if(expsum>expvalue)
		expsum=expsum+std::log(1.0+exp(expvalue-expsum));
	else
		expsum=expvalue+std::log(1.0+exp(expsum-expvalue));
}

inline double ITS_VES::exp_minus(double exp1,double exp2,bool& sign)
{
	if(exp1>exp2)
	{
		sign=false;
		return exp1+std::log(1.0-exp(exp2-exp1));
	}
	else
	{
		sign=true;
		return exp2+std::log(1.0-exp(exp1-exp2));
	}
}

inline double ITS_VES::exp_calc(double exp1,double exp2,bool sign1,bool sign2,bool& sign_fin)
{
	if(sign1==sign2)
	{
		sign_fin=sign1;
		return exp_add(exp1,exp2);
	}
	else if(sign1)
		return exp_minus(exp2,exp1,sign_fin);
	else
		return exp_minus(exp1,exp2,sign_fin);
}

inline std::vector<double> ITS_VES::omega2_alpha(const std::vector<std::vector<double> >& hessian,const std::vector<double>& dalpha)
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

inline std::vector<double> ITS_VES::omega2_alpha(const std::vector<double>& hessian,const std::vector<double>& dalpha)
{
	std::vector<double> result;
	for(unsigned i=0;i!=dalpha.size();++i)
		result.push_back(hessian[i]*dalpha[i]);
	return result;
}

void ITS_VES::output_fb()
{
	if(mcycle%fb_stride==0)
	{
		ofb.rewind();
		ofb.addConstantField("ITERATE_METHOD");
		if(is_ves)
			ofb.printField("ITERATE_METHOD","Variational");
		else
			ofb.printField("ITERATE_METHOD","Traditional");
		ofb.addConstantField("ITERATE_STEP").printField("ITERATE_STEP",int(mcycle));
		ofb.addConstantField("BOLTZMANN_CONSTANT").printField("BOLTZMANN_CONSTANT",kB);
		ofb.addConstantField("PESHIFT").printField("PESHIFT",peshift);
		for(unsigned i=0;i!=ntemp;++i)
		{
			string id;
			Tools::convert(i,id);
			ofb.printField("index",id);
			ofb.printField("temperature",int_temps[i]);
			ofb.printField("fb_value",fb[i]);
			if(is_ves)
				ofb.printField("iterate_fb",fb_iter[i]);
			else
				ofb.printField("norm_value",norml[i]);
			ofb.printField();
		}
		ofb.flush();
	}
	if(mcycle%fbtrj_stride==0)
	{
		ofbtrj.fmtField(" %f");
		ofbtrj.printField("step",int(mcycle));
		onormtrj.fmtField(" %f");
		onormtrj.printField("step",int(mcycle));
		opstrj.fmtField(" %f");
		opstrj.printField("step",int(mcycle));
		opstrj.printField("PESHIFT",peshift);
		for(unsigned i=0;i!=ntemp;++i)
		{
			string id;
			Tools::convert(i,id);
			string fbid="FB"+id;
			ofbtrj.printField(fbid,fb[i]);
			if(is_ves)
				onormtrj.printField(fbid,fb_iter[i]);
			else
				onormtrj.printField(fbid,norml[i]);
		}
		ofbtrj.printField();
		ofbtrj.flush();
		onormtrj.printField();
		onormtrj.flush();
		opstrj.printField();
		opstrj.flush();
	}
	if(bias_output&&mcycle%bias_stride==0)
		output_bias();
}

void ITS_VES::mw_merge_rbfb()
{
	multi_sim_comm.Sum(norm_step);
	if(comm.Get_rank()==0)
	{
		const unsigned& nw=multi_sim_comm.Get_size();
		std::vector<double> all_min_ener(nw,0);
		std::vector<double> all_rbfb(nw*rbfb.size(),0);
		multi_sim_comm.Allgather(rbfb,all_rbfb);
		multi_sim_comm.Allgather(min_ener,all_min_ener);

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
		double nave=std::log(double(nw));
		for(unsigned j=0;j!=rbfb.size();++j)
			rbfb[j]-=nave;
	}
	comm.Bcast(rbfb,0);
}

void ITS_VES::output_bias()
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
	for(unsigned i=1;i!=ntemp;++i)
		exp_added(sfb,fb[i]);

	for(unsigned i=0;i!=bias_bins;++i)
	{
		double energy=bias_min+d_pot*i;
		double senergy=energy+peshift;

		double _gfsum=-betak[0]*senergy+fb[0];
		double _bgfsum=_gfsum+std::log(betak[0]);
		for(unsigned j=1;j!=ntemp;++j)
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
