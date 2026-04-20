/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2024
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/*
 * Bekker-Wong soil plasticity constitutive law for MBDyn.
 *
 * Models terrain behavior during spacecraft/capsule landing using the
 * modified Bekker-Wong pressure-sinkage law with elastic unloading.
 * Suitable for granular soils: lunar/Martian regolith, terrestrial sand.
 *
 * Loading branch (z >= z_max, virgin loading or deeper penetration):
 *   F = A * (kc/b + kphi) * z^n
 *
 * Unloading/reloading branch (z < z_max):
 *   F = A * Ku * (z - z0)
 *   Ku = n * (kc/b + kphi) * z_max^(n-1)   [secant stiffness at z_max]
 *   z0 = z_max * (1 - 1/n)                  [permanent sinkage offset]
 *
 * Viscous damping (stabilizes impact transient):
 *   F_total = F_Bekker + c_damp * dz/dt
 *
 * Intended usage: pair with a Rod joint between the capsule leg-tip node
 * and a fixed ground node. One instance per leg — each tracks its own
 * permanent sinkage independently.
 *
 * Reference parameters (Bekker 1960, Wong 2008):
 *   Lunar regolith (loose): kc=1400, kphi=820000, n=1.0
 *   Lunar regolith (hard):  kc=5300, kphi=1740000, n=1.0
 *   Mars regolith:          kc=800,  kphi=140000,  n=0.9
 *   Terrestrial dry sand:   kc=990,  kphi=1528000, n=1.1
 *
 * Input syntax:
 *   constitutive law: <label>, 1D, bekker soil,
 *       [ , sign, { positive | negative | <real> }, ]
 *       kc, <kc>,           # cohesive modulus [N/m^(n+1)]
 *       kphi, <kphi>,       # frictional modulus [N/m^(n+2)]
 *       n, <n>,             # sinkage exponent [-], >= 0.5
 *       pad radius, <b>,    # effective contact radius [m]
 *       [ , pad area, <A>, ]  # override area (default: pi*b^2)
 *       damping, <c_damp>,  # viscous damping [N*s/m], >= 0
 *       [ , initial zmax, <z_max0>, ]   # restore plastic state on restart
 *       [ , prestrain, (DriveCaller)<prestrain> ]
 *       ;
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>
#include <algorithm>

#include "dataman.h"
#include "constltp_impl.h"

static const doublereal Z_FLOOR = 1.e-9;   /* sinkage floor [m] — avoids pow(0, n-1) */
static const doublereal M_PI_LOCAL = 3.14159265358979323846;

class BekkerSoilCL
: public ElasticConstitutiveLaw<doublereal, doublereal> {
private:
	/* Bekker-Wong soil parameters */
	doublereal m_dKc;       /* cohesive modulus [N/m^(n+1)]  */
	doublereal m_dKphi;     /* frictional modulus [N/m^(n+2)] */
	doublereal m_dN;        /* sinkage exponent [-]           */
	doublereal m_dB;        /* effective pad radius [m]       */
	doublereal m_dA;        /* effective contact area [m^2]   */
	doublereal m_dCdamp;    /* viscous damping [N*s/m]        */

	/* Sign convention: z = m_dSign * Epsilon (z > 0 means penetrating) */
	doublereal m_dSign;

	/* Plastic state — committed only in AfterConvergence() */
	doublereal m_dZmax;       /* maximum historical sinkage (committed) [m] */
	doublereal m_dZmax_trial; /* trial z_max accumulated during Newton-Raphson */
	doublereal m_dKu;         /* unloading stiffness [N/m] cached from z_max */
	doublereal m_dZ0;         /* permanent sinkage offset [m] cached from z_max */

	/* Contact state toggling (same pattern as HuntCrossleyCL) */
	bool m_bActive;
	bool m_bToggling;

	/* Precompute Bekker coefficient: K_bek = kc/b + kphi */
	doublereal KBek(void) const {
		return m_dKc / m_dB + m_dKphi;
	}

	/* Update cached unloading parameters from a given z_max value */
	void UpdateUnloadingCache(doublereal z_max) {
		if (z_max > Z_FLOOR) {
			m_dKu = m_dN * KBek() * std::pow(z_max, m_dN - 1.0);
			m_dZ0 = z_max - (m_dA * KBek() * std::pow(z_max, m_dN))
			              / (m_dA * m_dKu);
			/* Simplification: z0 = z_max*(1 - 1/n) for pure power law.
			 * The formula above is equivalent but numerically stable for
			 * any n >= 0.5. */
		} else {
			/* z_max too small: treat as elastic regime with tiny stiffness */
			m_dKu = m_dN * KBek() * std::pow(Z_FLOOR, m_dN - 1.0);
			m_dZ0 = 0.0;
		}
	}

public:
	BekkerSoilCL(const TplDriveCaller<doublereal> *pTplDC,
		doublereal dSign,
		doublereal dKc, doublereal dKphi, doublereal dN,
		doublereal dB, doublereal dA,
		doublereal dCdamp,
		doublereal dZmax0)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pTplDC, 0.),
	m_dKc(dKc), m_dKphi(dKphi), m_dN(dN),
	m_dB(dB), m_dA(dA), m_dCdamp(dCdamp),
	m_dSign(dSign),
	m_dZmax(dZmax0), m_dZmax_trial(dZmax0),
	m_dKu(0.), m_dZ0(0.),
	m_bActive(false), m_bToggling(false)
	{
		UpdateUnloadingCache(m_dZmax);
	};

	virtual ~BekkerSoilCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		/* VISCOELASTIC because Update() uses both Eps and EpsPrime */
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		SAFENEWWITHCONSTRUCTOR(pCL, BekkerSoilCL,
			BekkerSoilCL(pGetDriveCaller()->pCopy(),
				m_dSign,
				m_dKc, m_dKphi, m_dN,
				m_dB, m_dA,
				m_dCdamp,
				m_dZmax));  /* copy current plastic state */
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "bekker soil"
			<< ", sign, "        << m_dSign
			<< ", kc, "          << m_dKc
			<< ", kphi, "        << m_dKphi
			<< ", n, "           << m_dN
			<< ", pad radius, "  << m_dB
			<< ", pad area, "    << m_dA
			<< ", damping, "     << m_dCdamp
			<< ", initial zmax, "<< m_dZmax  /* preserve plastic state on restart */
			<< ", ";
		ElasticConstitutiveLaw<doublereal, doublereal>::Restart_int(out);
		return out;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		/* Strain relative to prestrain drive */
		ConstitutiveLaw<doublereal, doublereal>::Epsilon =
			Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		/* Sinkage: z > 0 means the leg is penetrating the soil */
		doublereal z  = m_dSign * ConstitutiveLaw<doublereal, doublereal>::Epsilon;
		doublereal zp = m_dSign * ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime;

		/* ---- No contact: leg above ground ---- */
		if (z <= 0.) {
			if (m_bActive && !m_bToggling) {
				m_bToggling = true;   /* will deactivate at AfterConvergence */
			}
			ConstitutiveLaw<doublereal, doublereal>::F        = 0.;
			ConstitutiveLaw<doublereal, doublereal>::FDE      = 0.;
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = 0.;
			return;
		}

		/* ---- Contact active ---- */
		if (!m_bActive && !m_bToggling) {
			m_bToggling = true;       /* will activate at AfterConvergence */
		}

		const doublereal K_bek = KBek();

		/* Update trial z_max (not committed — only updated at AfterConvergence) */
		if (z > m_dZmax_trial) {
			m_dZmax_trial = z;
		}

		doublereal F_el, dFdz;

		if (z >= m_dZmax) {
			/* ===== Loading branch: virgin or deeper penetration ===== *
			 * F = A * K_bek * z^n
			 * dF/dz = A * K_bek * n * z^(n-1)                         */
			doublereal z_safe = std::max(z, Z_FLOOR);
			doublereal zn     = std::pow(z_safe, m_dN);
			doublereal znm1   = std::pow(z_safe, m_dN - 1.0);
			F_el  = m_dA * K_bek * zn;
			dFdz  = m_dA * K_bek * m_dN * znm1;

		} else {
			/* ===== Unloading/reloading branch ===== *
			 * F = A * Ku * (z - z0)  [linear elastic from committed z_max] *
			 * Ku and z0 are cached from the last committed z_max.           */
			doublereal dz = z - m_dZ0;
			if (dz < 0.) {
				dz = 0.;   /* cannot pull soil upward (no adhesion) */
			}
			F_el  = m_dA * m_dKu * dz;
			dFdz  = m_dA * m_dKu;
		}

		/* Viscous damping (stabilizes impact transient) */
		doublereal F_damp = m_dCdamp * zp;

		/* Total force and Jacobians.
		 *
		 * FDE (tangent stiffness dF/dEps):
		 *   F = m_dSign * (F_el(z) + F_damp)  with  z = m_dSign * Eps
		 *   dF/dEps = m_dSign * dFel/dz * dz/dEps = m_dSign * dFdz * m_dSign = dFdz
		 * So FDE = dFdz regardless of sign, always >= 0. */
		ConstitutiveLaw<doublereal, doublereal>::F        = m_dSign * (F_el + F_damp);
		ConstitutiveLaw<doublereal, doublereal>::FDE      = dFdz;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = m_dSign * m_dCdamp;
	};

	virtual void AfterConvergence(const doublereal& Eps,
		const doublereal& EpsPrime = 0.)
	{
		/* ---- Commit plastic state ---- *
		 * If the trial z_max exceeds the committed value, the soil has
		 * permanently deformed further. Update cached unloading parameters. */
		if (m_dZmax_trial > m_dZmax) {
			m_dZmax = m_dZmax_trial;
			UpdateUnloadingCache(m_dZmax);
		}
		/* Reset trial accumulator for the next timestep */
		m_dZmax_trial = m_dZmax;

		/* ---- Commit contact state toggle ---- */
		if (m_bToggling) {
			if (m_bActive) {
				m_bActive   = false;
				/* Leg lifted off: keep m_dZmax (crater persists).
				 * Reset trial to committed value so next contact
				 * starts fresh from the existing crater. */
				m_dZmax_trial = m_dZmax;
			} else {
				m_bActive = true;
			}
			m_bToggling = false;
		}
	};
};

/* ---- Reader ---- */
struct BekkerSoilCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cerr(
				"BekkerSoilCL:\n"
				"        bekker soil,\n"
				"                [ , sign, { negative | positive | <sign> } , ]\n"
				"                kc,         <kc>,       # cohesive modulus [N/m^(n+1)]\n"
				"                kphi,       <kphi>,     # frictional modulus [N/m^(n+2)]\n"
				"                n,          <n>,        # sinkage exponent [-], >= 0.5\n"
				"                pad radius, <b>,        # effective contact radius [m]\n"
				"                [ , pad area, <A>, ]    # override area (default: pi*b^2)\n"
				"                damping,    <c_damp>,   # viscous damping [N*s/m], >= 0\n"
				"                [ , initial zmax, <zmax0>, ]  # restore plastic state\n"
				"                [ , prestrain, (DriveCaller)<prestrain> ]\n"
				"\n"
				"  Reference Bekker parameters (kc [N/m^(n+1)], kphi [N/m^(n+2)], n):\n"
				"    Lunar regolith (loose): kc=1400,  kphi=820000,  n=1.0\n"
				"    Lunar regolith (hard):  kc=5300,  kphi=1740000, n=1.0\n"
				"    Mars regolith:          kc=800,   kphi=140000,  n=0.9\n"
				"    Terrestrial dry sand:   kc=990,   kphi=1528000, n=1.1\n"
				<< std::endl);

			if (!HP.IsArg()) {
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}
		}

		/* sign (default: negative — compression is negative Epsilon) */
		doublereal dSign = -1.;
		if (HP.IsKeyWord("sign")) {
			if (HP.IsKeyWord("positive")) {
				dSign = 1.;
			} else if (HP.IsKeyWord("negative")) {
				dSign = -1.;
			} else {
				doublereal d = HP.GetReal();
				if (d == 0.) {
					silent_cerr("BekkerSoilCLR: invalid sign " << d
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				dSign = copysign(1., d);
			}
		}

		/* kc */
		if (!HP.IsKeyWord("kc")) {
			silent_cerr("BekkerSoilCLR: \"kc\" expected"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dKc = HP.GetReal();
		if (dKc < 0.) {
			silent_cerr("BekkerSoilCLR: invalid \"kc\" " << dKc
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* kphi */
		if (!HP.IsKeyWord("kphi")) {
			silent_cerr("BekkerSoilCLR: \"kphi\" expected"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dKphi = HP.GetReal();
		if (dKphi <= 0.) {
			silent_cerr("BekkerSoilCLR: invalid \"kphi\" " << dKphi
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* n */
		if (!HP.IsKeyWord("n")) {
			silent_cerr("BekkerSoilCLR: \"n\" expected"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dN = HP.GetReal();
		if (dN < 0.5) {
			silent_cerr("BekkerSoilCLR: invalid \"n\" " << dN
				<< " (must be >= 0.5)"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* pad radius */
		if (!HP.IsKeyWord("pad" "radius")) {
			silent_cerr("BekkerSoilCLR: \"pad radius\" expected"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dB = HP.GetReal();
		if (dB <= 0.) {
			silent_cerr("BekkerSoilCLR: invalid \"pad radius\" " << dB
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* pad area (optional — defaults to pi*b^2) */
		doublereal dA = M_PI_LOCAL * dB * dB;
		if (HP.IsKeyWord("pad" "area")) {
			dA = HP.GetReal();
			if (dA <= 0.) {
				silent_cerr("BekkerSoilCLR: invalid \"pad area\" " << dA
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		/* damping */
		if (!HP.IsKeyWord("damping")) {
			silent_cerr("BekkerSoilCLR: \"damping\" expected"
				" at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal dCdamp = HP.GetReal();
		if (dCdamp < 0.) {
			silent_cerr("BekkerSoilCLR: invalid \"damping\" " << dCdamp
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* initial zmax (optional — for restarting with existing plastic state) */
		doublereal dZmax0 = 0.;
		if (HP.IsKeyWord("initial" "zmax")) {
			dZmax0 = HP.GetReal();
			if (dZmax0 < 0.) {
				silent_cerr("BekkerSoilCLR: invalid \"initial zmax\" " << dZmax0
					<< " at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		/* Prestrain drive (optional) */
		TplDriveCaller<doublereal> *pTplDC = GetPreStrain<doublereal>(pDM, HP);

		SAFENEWWITHCONSTRUCTOR(pCL, BekkerSoilCL,
			BekkerSoilCL(pTplDC, dSign,
				dKc, dKphi, dN,
				dB, dA,
				dCdamp,
				dZmax0));

		return pCL;
	};
};

/* ---- Module registration ---- */
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new BekkerSoilCLR;
	if (!SetCL1D("bekker" "soil", rf1D)) {
		delete rf1D;
		silent_cerr("BekkerSoilCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}
	return 0;
}
