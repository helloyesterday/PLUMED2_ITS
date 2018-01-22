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
	double cv_energy;
	double shift_energy;
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
	std::vector<double> norml;

	std::vector<double> gU;
	std::vector<double> gf;
	std::vector<double> bgf;
	std::vector<double> rbfb;
	std::vector<double> pot_aver;
	std::vector<double> pot_norm;

	std::vector<double> energy_mean;
	std::vector<double> energy_dev;
	std::vector<double> fb_iter;

	std::vector<double> rw_temp;
	std::vector<double> rw_beta;
	std::vector<double> rw_factor;
	std::vector<double> rw_dth;
	std::vector<double> rw_dtl;
	std::vector<unsigned> rw_rctid;
	std::vector<unsigned> pot_dis;
	std::vector<std::string> rw_tstr;

	std::vector<double> energy_record;
	std::vector<double> gfsum_record;
	std::vector<double> bgfsum_record;
	std::vector<double> effpot_record;
	std::vector<double> bias_record;
	std::vector<double> force_record;
	std::vector<double> bias_ratio;

	// The average of 1st order drivative of bias potential in target distribution
	std::vector<double> deriv_ps;
	// The average of 2nd order drivative of bias potential in target distribution
	std::vector<double> deriv2_ps;

	unsigned long step;
	unsigned norm_step;
	unsigned mcycle;
	unsigned iter_limit;
	double fb_init;
	double fb_bias;
	double rb_fac1;
	double rb_fac2;
	double fb0;
	double bias_min;
	double bias_max;
	double d_pot;
	double step_size;

	bool equiv_temp;
	bool is_const;
	bool rw_output;
	bool is_ves;
	bool read_norm;
	bool only_1st;
	bool use_mw;
	bool bias_output;
	bool rbfb_output;
	bool is_debug;
	bool potdis_output;
	bool bias_linked;
	bool only_bias;
	bool use_fixed_peshift;
	bool is_unlinear;
	bool is_set_temps;
	bool is_set_ratios;
	bool is_norm_rescale;
	bool is_direct;
	bool read_fb;
	bool read_iter;
	bool fbtrj_output;
	bool is_rw_rct;
	//~ bool norm_output;
	//~ bool peshift_output;
	unsigned update_step;
	unsigned rw_stride;
	unsigned rw_start;
	unsigned start_cycle;
	unsigned fb_stride;
	unsigned fbtrj_stride;
	unsigned bias_stride;
	unsigned bias_bins;
	unsigned potdis_num;
	unsigned potdis_step;
	unsigned potdis_update;
	unsigned nbiases;
	unsigned rctid;

	std::string fb_file;
	std::string fb_trj;
	std::string norm_trj;
	std::string iter_trj;
	std::string deriv_trj;
	std::string peshift_trj;
	std::string rw_file;
	std::string fb_input;
	std::string target_file;
	std::string bias_file;
	std::string debug_file;
	std::string potdis_file;
	std::string rbfb_file;

	OFile ofb;
	OFile ofbtrj;
	//~ OFile onormtrj;
	//~ OFile oderiv;
	//~ OFile opstrj;
	OFile ofw;
	OFile obias;
	OFile orbfb;
	OFile odebug;
	OFile opotdis;

	double gfsum;
	double bgfsum;
	double min_ener;
	double pot_min;
	double pot_max;
	double pot_bin;
	double tot_bias;
	double sim_dtl;
	double sim_dth;

	unsigned ntarget;
	double ener_min;
	double ener_max;
	double dU;
	double dvp2_complete;
	
	std::vector<Bias*> bias_pntrs_;

	Value* valueEnergy;
	Value* valueEff;
	Value* valueForce;
	Value* valueRct;
	Value* valueRBias;
	Value* valueRwfb;

	inline void update_rbfb();
	inline void update_rbfb_direct();
	
	void fb_iteration();
	void fb_variational();
	void output_fb();
	void mw_merge_rbfb();
	void output_bias();
	void setupOFile(std::string& file_name, OFile& ofile, const bool multi_sim_single_files);
	unsigned read_fb_file(const std::string& fname,double& _kB,double& _peshift);
	double calc_deriv2();
	unsigned find_rw_id(double rwtemp,double& dtl,double& dth);
	double find_rw_fb(unsigned rwid,double dtl,double dth);
	double find_rw_fb(double rwtemp);

	inline void coe_rescale(double shift,std::vector<double>& coe);
	inline void fb_rescale(double shift){coe_rescale(shift,fb);}
	inline void iter_rescale(double shift){coe_rescale(shift,fb_iter);}
	void change_peshift(double new_shift);

	std::vector<double> omega2_alpha(
		const std::vector<std::vector<double> >& hessian,
		const std::vector<double>& dalpha);
	std::vector<double> omega2_alpha(
		const std::vector<double>& hessian,
		const std::vector<double>& dalpha);

public:
	explicit ITS_Bias(const ActionOptions&);
	~ITS_Bias();
	void calculate();
	double calc_bias(double _beta){return -gfsum/_beta-shift_energy;}
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
