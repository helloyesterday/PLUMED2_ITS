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
#ifndef __PLUMED_bias_ITS_Bias_h
#define __PLUMED_bias_ITS_Bias_h

#include "Bias.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include "tools/Tools.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <ctime>
#include <cmath>

namespace PLMD{
namespace bias{

class ITS_Bias: public Bias
{
private:
	const unsigned update_start;
	unsigned nreplica;
	double templ;
	double temph;
	double ratiol;
	double ratioh;
	double kB;
	double kBT;
	double kBT_target;
	double beta0;
	double beta_target;
	double target_temp;
	double sim_temp;
	double peshift;
	double energy;
	double shift_pot;
	double cv_energy;
	double input_energy;
	double eff_energy;
	double eff_factor;
	double bias_energy;
	double bias_force;
	double rct;
	double rbias;

	std::vector<double> int_temps;
	std::vector<double> int_ratios;
	std::vector<double> betak;

	std::vector<double> fb;
	std::vector<double> fb_rct;
	std::vector<double> norml;
	std::vector<double> fb_ratios;
	std::vector<double> rbfb_ratios;

	std::vector<double> gU;
	std::vector<double> gE;
	std::vector<double> gf;
	std::vector<double> bgf;
	std::vector<double> rbfb;
	//~ std::vector<double> rbfb2;
	std::vector<double> rbzb;
	std::vector<double> pot_aver;
	std::vector<double> pot_norm;
	std::vector<double> zbe;
	std::vector<double> zbu;

	std::vector<double> energy_mean;
	std::vector<double> energy_dev;
	std::vector<double> fb_iter;

	std::vector<double> rw_temp;
	std::vector<double> rw_beta;
	std::vector<double> rw_factor;
	std::vector<double> rw_dth;
	std::vector<double> rw_dtl;
	std::vector<double> rw_rct;
	std::vector<double> rw_fb_ratios;
	std::vector<double> rw_rctid;
	std::vector<unsigned> pot_dis;
	std::vector<unsigned> gmin_id;
	std::vector<unsigned> gmax_id;
	std::vector<std::string> rw_tstr;

	std::vector<double> energy_record;
	std::vector<double> gfsum_record;
	std::vector<double> bgfsum_record;
	std::vector<double> effpot_record;
	std::vector<double> bias_record;
	std::vector<double> force_record;
	std::vector<double> bias_ratio;
	std::vector<double> peshift_ratio;
	std::vector<double> peavg;
	std::vector<double> peavg2;
	std::vector<double> pedev;
	std::vector<double> penorm;
	std::vector<double> gauss_mean;
	std::vector<double> gauss_sd;
	std::vector<double> gauss_sd2;
	std::vector<double> gauss_weight0;
	std::vector<double> gauss_weight1;
	std::vector<double> energy_range;
	std::vector<double> eff_p0;
	std::vector<double> eff_p1;
	std::vector<std::vector<double> > gauss_values;

	unsigned long step;
	unsigned norm_step;
	unsigned mcycle;
	unsigned iter_limit;
	unsigned ema_days;
	//~ unsigned pre_iter_steps;
	unsigned avg_step;
	double fb_init;
	double fb_bias;
	double rb_fac1;
	double rb_fac2;
	double fb0;
	double bias_min;
	double bias_max;
	double d_pot;
	double temp_ratio_energy;
	double anneal_rate;
	double kl_rate;
	double int_dE;
	//~ double neighbor_ratio;

	bool equiv_temp;
	bool is_const;
	bool rw_output;
	bool read_norm;
	bool only_1st;
	bool use_mw;
	bool bias_output;
	bool rbfb_output;
	bool zke_output;
	bool zku_output;
	bool is_debug;
	bool potdis_output;
	bool bias_linked;
	bool only_bias;
	bool is_read_ratio;
	bool auto_peshift;
	bool is_unlinear;
	bool is_set_temps;
	bool is_set_ratios;
	bool is_norm_rescale;
	bool read_fb;
	bool read_iter;
	bool fbtrj_output;
	bool rct_output;
	bool no_bias_rct;
	bool read_zbu;
	bool rw_partition;
	bool do_anneal;
	bool use_ema;
	bool abs_ratio;
	bool output_pegauss;
	bool do_calc_pegauss;
	//~ bool norm_output;
	//~ bool peshift_output;
	
	unsigned update_step;
	unsigned start_cycle;
	unsigned fb_stride;
	unsigned fbtrj_stride;
	unsigned rct_stride;
	unsigned bias_stride;
	unsigned bias_bins;
	unsigned potdis_num;
	unsigned potdis_step;
	unsigned potdis_update;
	unsigned nbiases;
	unsigned int_num;
	unsigned dE_num;
	unsigned pegauss_update;

	std::string fb_file;
	std::string fb_trj;
	std::string norm_trj;
	std::string iter_trj;
	std::string deriv_trj;
	std::string peshift_trj;
	std::string fb_input;
	std::string target_file;
	std::string bias_file;
	std::string debug_file;
	std::string potdis_file;
	std::string rbfb_file;
	std::string zke_file;
	std::string zku_file;
	std::string rct_file;
	std::string pegauss_file;

	OFile ofb;
	OFile ofbtrj;
	OFile ofw;
	OFile obias;
	OFile orbfb;
	OFile ozke;
	OFile ozku;
	OFile orct;
	OFile odebug;
	OFile opotdis;
	OFile opegauss;

	double gfsum;
	double bgfsum;
	double min_ener;
	double pot_min;
	double pot_max;
	double pot_bin;
	double tot_bias;
	//~ double rbfb_ratio;
	double logN;
	double Nlogr;
	double ratio_norm;
	double fb_ratio0;
	double rctid;
	double wemat;
	double wemay;

	unsigned ntarget;
	double ener_min;
	double ener_max;
	double dU;
	double avg_shift;
	double gauss_limit;
	
	std::vector<Bias*> bias_pntrs_;
	
	//~ Value* valueEnergy;
	Value* valuePot;
	Value* valueEff;
	Value* valueForce;
	Value* valueRBias;
	Value* valueRwfb;
	std::vector<Value*> valueRwbias;

	inline void update_rbfb();
	inline void update_peavg();
	inline void update_rw_rct(bool use_new_ratios);
	inline void update_ratios(double new_energ);
	inline bool kl_divergence(const std::vector<double>& p0,const std::vector<double>& p1);
	bool anneal_judge(double new_energy);

	void calc_pe_vag();
	void fb_iteration();
	void zbu_iteration();
	void fb_partition();
	void output_fb();
	void mw_merge_data();
	void output_bias();
	void calc_gaussian();
	void update_gaussian_grid();
	void annealing_fb();
	void setupOFile(std::string& file_name, OFile& ofile, const bool multi_sim_single_files);
	void setupOFiles(std::vector<std::string>& fnames, std::vector<OFile*>& OFiles, const bool multi_sim_single_files);
	unsigned read_fb_file(const std::string& fname,double& _kB,double& _peshift);
	double find_rw_id(double rwtemp);
	double find_rw_fb(double real_id);

	inline void coe_rescale(double shift,std::vector<double>& coe);
	inline void fb_rescale(double shift){coe_rescale(shift,fb);}
	inline void iter_rescale(double shift){coe_rescale(shift,fb_iter);}
	void change_peshift(double new_shift);
	void set_peshift_ratio();

public:
	explicit ITS_Bias(const ActionOptions&);
	~ITS_Bias();
	void calculate();
	double calc_bias(double _beta) const
		{return -gfsum/_beta-input_energy;}
	double calc_rct(double _beta,double _fb,double _fb_ratio) const
		{return (_fb_ratio-_fb)/_beta;}
	//~ double calc_rct(double _beta,double _fb) const
		//~ {return (-logN-_fb)/_beta;}
	static void registerKeywords(Keywords& keys);
};

inline void ITS_Bias::coe_rescale(double shift,std::vector<double>& coe)
{
	double coeshift=coe[0]+betak[0]*shift;
	for(unsigned i=0;i!=coe.size();++i)
		coe[i]+=betak[i]*shift-coeshift;
}

}
}

#endif
