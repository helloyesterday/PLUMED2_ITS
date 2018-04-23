/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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


namespace PLMD {
namespace bias {

Bias::Bias(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  outputForces(getNumberOfArguments(),0.0),
  ex_bias_linked_(false),
  is_set_rct(false),
  ex_bias_pntr_(NULL),
  extra_bias_ratio_(0)
{
  addComponentWithDerivatives("bias");
  componentIsNotPeriodic("bias");
  valueBias=getPntrToComponent("bias");

  if(getStride()>1) {
    log<<"  multiple time step "<<getStride()<<" ";
    log<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  }

  turnOnDerivatives();
}

void Bias::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
}

void Bias::apply() {
  const unsigned noa=getNumberOfArguments();
  const unsigned ncp=getNumberOfComponents();
  
  // apply the ratio of extra bias
  if(ex_bias_linked_)
  {
    for(unsigned i=0; i<noa; ++i)
    {
      double opf=outputForces[i];
      outputForces[i]*=1.0+extra_bias_ratio_;
      valueBias->addDerivative(i,-opf*extra_bias_ratio_);
	}
    extra_bias_ratio_=0;
  }

  if(onStep()) {
    double gstr = static_cast<double>(getStride());
    for(unsigned i=0; i<noa; ++i) {
      getPntrToArgument(i)->addForce(gstr*outputForces[i]);
    }
  }

  // additional forces on the bias component
  std::vector<double> f(noa,0.0);
  std::vector<double> forces(noa);

  bool at_least_one_forced=false;
  for(unsigned i=0; i<ncp; ++i) {
    if(getPntrToComponent(i)->applyForce(forces)) {
      at_least_one_forced=true;
      for(unsigned j=0; j<noa; j++) f[j]+=forces[j];
    }
  }

  if(at_least_one_forced && !onStep()) error("you are biasing a bias with an inconsistent STRIDE");

  if(at_least_one_forced) for(unsigned i=0; i<noa; ++i) {
      getPntrToArgument(i)->addForce(f[i]);
    }

}

void Bias::setRctComponent(const std::string& rct_name){
  addComponent(rct_name);
  componentIsNotPeriodic(rct_name);
  valueRct=getPntrToComponent(rct_name);
  valueRct->set(0);
  is_set_rct=true;
}

void Bias::linkExternalBias(Bias* ex_bias_pntr_in)
{
  if(ex_bias_pntr_==NULL){
    ex_bias_pntr_ = ex_bias_pntr_in;
  }
  else {
    std::string err_msg = "bias " + getLabel() + " of type " + getName() +
      " has already been linked with ITS pointor " + ex_bias_pntr_->getLabel() +
      " of type " + ex_bias_pntr_->getName() + ". You cannot link two ITS pointer to the same bias.";
    plumed_merror(err_msg);
  }
  ex_bias_linked_=true;
}

}
}


