#include "nec_results.h"
#include "nec_debug.h"

void nec_antenna_input::set_input(int tag, int segment, nec_complex voltage, nec_complex current, nec_complex impedance, nec_complex admittance, nec_float power)
{
  DEBUG_TRACE("nec_antenna_input::set_input(" << tag << "," << segment << "," << voltage << "," << current  << ",");
  DEBUG_TRACE("                             " << impedance << "," << admittance << "," << power << ")");
  n_items++;
  _tag.push_back(tag);
  _segment.push_back(segment);
  _voltage.push_back(voltage);
  _current.push_back(current);
  _impedance.push_back(impedance);
  _admittance.push_back(admittance);
  _power.push_back(power);
}

vector<nec_complex>& nec_antenna_input::get_impedance()
{
  DEBUG_TRACE("get_impedance()");
  DEBUG_TRACE("element count :" << _impedance.size());
 return _impedance;
}
                


void nec_near_field_pattern::set_input(nec_float x, nec_float y, nec_float z, nec_complex field_x, nec_complex field_y, nec_complex field_z)
{
  DEBUG_TRACE("nec_near_field_pattern::set_input(" << x << "," << y << "," << z << "," << field_x  << ")");
  n_items++;
  _x.push_back(x);
  _y.push_back(y);
  _z.push_back(z);
  _field_x.push_back(field_x);
  _field_y.push_back(field_y);
  _field_z.push_back(field_z);
}                
