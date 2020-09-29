#include "openmc/tallies/filter_legendre.h"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
LegendreFilter::from_xml(pugi::xml_node node)
{
  this->set_order(std::stoi(get_node_value(node, "order")));

  auto angle = get_node_value(node, "angle");
  if (angle == "mu") {
    this->set_angle(LegendreAngle::mu);
  } else if (angle == "azimuthal") {
    this->set_angle(LegendreAngle::azimuthal);
  } else if (angle == "polar") {
    this->set_angle(LegendreAngle::polar);
  } else {
    throw std::runtime_error{"Angle for LegendreFilter must be 'mu', "
                             "'azimuthal', or 'polar'"};
  }
}

void
LegendreFilter::set_order(int order)
{
  if (order < 0) {
    throw std::invalid_argument{"Legendre order must be non-negative."};
  }
  order_ = order;
  n_bins_ = order_ + 1;
}

void
LegendreFilter::set_angle(LegendreAngle angle)
{
  angle_ = angle;
}

void
LegendreFilter::get_all_bins(const Particle& p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  // Get the angle of interest.
  double p_angle;
  if (angle_ == LegendreAngle::mu) {
    p_angle = p.mu_;
  } else if (angle_ == LegendreAngle::azimuthal) {
    if (abs(p.u()[0]) > 1e-12)
      p_angle = atan(p.u()[1] / p.u()[0]) / PI * 2 - 1;
    else
      p_angle = (p.u()[0] >= 0) ? PI / 2 : ((p.u()[0] < 0) ? -PI / 2 : 0);
  } else {
    p_angle = acos(p.u()[2]) / PI * 2 - 1;
  }

  std::vector<double> wgt(n_bins_);
  calc_pn_c(order_, p_angle, wgt.data());
  for (int i = 0; i < n_bins_; i++) {
    match.bins_.push_back(i);
    match.weights_.push_back(wgt[i]);
  }
}

void
LegendreFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);

  if (angle_ == LegendreAngle::mu) {
    write_dataset(filter_group, "angle", "mu");
  } else if (angle_ == LegendreAngle::azimuthal) {
    write_dataset(filter_group, "angle", "azimuthal");
  } else {
    write_dataset(filter_group, "angle", "polar");
  }
}

std::string
LegendreFilter::text_label(int bin) const
{
  if (angle_ == LegendreAngle::mu) {
    return fmt::format("Legendre expansion in mu, P{}", bin);
  } else if (angle_ == LegendreAngle::azimuthal) {
    return fmt::format("Legendre expansion in the azimuthal angle, P{}", bin);
  } else {
    return fmt::format("Legendre expansion in the polar angle, P{}", bin);
  }
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_legendre_filter_get_order(int32_t index, int* order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the order.
  *order = filt->order();
  return 0;
}

extern "C" int
openmc_legendre_filter_get_angle(int32_t index, int* angle)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the order.
  *angle = static_cast<int>(filt->angle());
  return 0;
}

extern "C" int
openmc_legendre_filter_set_order(int32_t index, int order)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_order(order);
  return 0;
}

extern "C" int
openmc_legendre_filter_set_angle(int32_t index, const int* angle)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<LegendreFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not a legendre filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_angle(static_cast<LegendreAngle>(*angle));
  return 0;
}

} // namespace openmc
