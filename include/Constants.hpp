#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>

namespace Constants {
	template <typename T>
	inline T pi(void){return boost::math::constants::pi<T>();}
	template <typename T>
	inline T root_pi(void) { return boost::math::constants::root_pi<T>(); }
	template <typename T>
	inline T half_pi(void) { return boost::math::constants::half_pi<T>(); }
	template <typename T>
	inline T pi_2(void) { return half_pi<T>(); }
	template <typename T>
	inline T pi_3(void){return boost::math::constants::pi<T>()/T(3);}
	template <typename T>
	inline T pi_4(void){return boost::math::constants::pi<T>()/T(4);}
	template <typename T>
	inline T twopi(void){return T(2)*boost::math::constants::pi<T>();}

	template <typename T>
	inline T Eh(void){return T(27.2113845);} // eV/Hartree
	template <typename T>
	inline T a0(void) { return T(0.5291772108); } // Ang/bohr
	template <typename T>
	inline T icm(void) { return T(8065.54445); } // icm/eV
	template <typename T>
	inline T hc(void) { return T(1239.84190604789); } // 197.326968*2*pi eV nm
	template <typename T>
	inline T C_cmPs(void) { return T(2.99792458e10); } //2.99792458e10 // cm/s
	template <typename T>
	inline T C_nmPfs(void) { return T(2.99792458e2); } //2.99792458e2 // nm/fs
	template <typename T>
	inline T MC2(void) { return T(931.494028e6); } // eV [u of atomic mass in eV]
	template <typename T>
	inline T mc2inau(void) { return MC2<T>()/Eh<T>(); }
	template <typename T>
	inline T amuPe_mass(void) { return T(5.4857990943e-4); } // electron mass in u of atomic mass
	template <typename T>
	inline T hbar(void) { return T(6.58211915e-16); }
	template <typename T>
	inline T kb(void) { return T(8.617343e-5); } // eV / K
	template <typename T>
	inline T fsPau(void) { return T(.02418884326505); } // fs / au
	template <typename T>
	inline T cinau(void) { return T(2997.92458)*fsPau<T>()/a0<T>(); } // C in units bohr/au
	template <typename T>
	inline T icmPau(void) { return icm<T>() * Eh<T>(); } // icm/hartree
	template <typename T>
	inline T amuPau(void) { return T(9.1093826/1.66053886)*T(1e-4); } // unified AMU / au_mass
	template <typename T>
	inline T e2P4pieps0(void) {return T(14.3996445); } //eV Ang
	template <typename T>
	inline T auPe2P4pieps0(void) { return (e2P4pieps0<T>()/Eh<T>()/a0<T>()); }  // hartree borh
	template <typename T>
	inline T aufor10PW(void) { return T(0.5336); } // Atomic unit of field for 10^16 W cm^-1
	template <typename T>
	inline T auenergy(void) { return T(183.631526); }// 4pi eps_0 * E_au^2 in units of  eV ang^-3
	template <typename T>
	inline T mbarn(void) {return T(1e-24); } // Mb/cm^2
	template <typename T>
	inline T auIntensity(void) { return T(3.55e16); }  // W/cm^2
	template <typename T>
	inline T aupolarizability(void) { return std::pow(a0<T>(),int(3)); } // ang^3
};







#endif

/*

Upinau(x,y) = 8*M_PI/137 * x/4/(y**2) // all in atomic units 
 = 0.1835 * x/ 4/y
 = 200 * x so for 3micron light, 1au of Up needs 0.005 au of intensity 
*/


