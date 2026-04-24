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
 * Surface impact element for MBDyn — Bekker-Wong plasticity + optional
 * hydrodynamic contribution.
 *
 * ---- Normal force (always active when node is in contact) ----
 *
 * 1) Bekker-Wong with elastic-plastic load/unload (MANDATORY):
 *
 *    Loading   (z >= z_max, virgin or deepening penetration):
 *      F_BEK = A * (kc/b + kphi) * z^n
 *
 *    Unloading / reloading (z < z_max):
 *      F_BEK = A * Ku * (z - z0)
 *      Ku  = n * (kc/b + kphi) * z_max^(n-1)   [secant stiffness at z_max]
 *      z0  = z_max * (1 - 1/n)                  [permanent sinkage offset]
 *
 *    z_max is the maximum historical penetration, committed per node at each
 *    converged timestep via AfterConvergence().
 *
 * ---- Total normal force ----
 *    FN = F_BEK   (clamped >= 0, no adhesion)
 *
 * ---- Coulomb friction (OPTIONAL, keyword "friction") ----
 *    F_fric = -mu_s * FN * (VT / |VT|)
 *    where VT is the tangential velocity component.
 *    Uses the TOTAL normal force FN.
 *
 * ---- Plasticity state management ----
 *    Each listed node has its own z_max (independent crater depth).
 *    Within each timestep, z_max_trial is updated monotonically in AssRes.
 *    AfterConvergence() commits z_max_trial to z_max and recomputes
 *    Ku and z0 for the next step.
 *
 * Input syntax:
 *   user defined: <label>, surface impact,
 *       plane normal, <nx>, <ny>, <nz>,   # outward unit normal (auto-normalised)
 *       plane point,  <px>, <py>, <pz>,   # any point on the plane
 *       kc,         <kc>,                 # Bekker cohesive modulus [N/m^(n+1)]
 *       kphi,       <kphi>,               # Bekker frictional modulus [N/m^(n+2)]
 *       n,          <n>,                  # sinkage exponent [-], >= 0.5
 *       pad radius, <b>,                  # smaller sinkage dimension [m]
 *       [ , area, <anor>, ]               # contact area [m^2] (default: pi*b^2)
 *       [ , friction, <mu_s>, ]           # Coulomb coefficient [-]
 *       nodes, <N>,                       # number of nodes
 *           <label_1>, ..., <label_N>     # structural node labels
 *       ;
 *
 * In control data:
 *   loadable elements: 1;
 *
 * Registration key: "surface impact"
 *
 * Reference Bekker parameters (kc [N/m^(n+1)], kphi [N/m^(n+2)], n):
 *   Lunar regolith (loose): kc=1400,  kphi=820000,  n=1.0
 *   Lunar regolith (hard):  kc=5300,  kphi=1740000, n=1.0
 *   Mars regolith:          kc=800,   kphi=140000,  n=0.9
 *   Terrestrial dry sand:   kc=990,   kphi=1528000, n=1.1
 */

#include "mbconfig.h"

#include <cmath>
#include <limits>
#include <vector>

#include "dataman.h"
#include "userelem.h"
#include "strnode.h"

static const doublereal Z_FLOOR = 1.e-9;
static const doublereal M_PI_L  = 3.14159265358979323846;

/* ========================================================================= */
class SurfaceImpactElem : public UserDefinedElem {
private:
	/* ---- plane ---- */
	Vec3 m_Normal;
	Vec3 m_PlanePoint;

	/* ---- Bekker-Wong (mandatory) ---- */
	doublereal m_dKc;        /* cohesive modulus [N/m^(n+1)] */
	doublereal m_dKphi;      /* frictional modulus [N/m^(n+2)] */
	doublereal m_dNbek;      /* sinkage exponent [-] */
	doublereal m_dB;         /* pad radius / smaller dimension [m] */
	doublereal m_dANOR;      /* contact area [m^2] */
	doublereal m_dKplastic;  /* plastic ratio: z0 = k_p * z_max, 0 < k_p < 1
	                          * fraction of z_max that remains as permanent sinkage.
	                          * k_p=0 → elastic (no dissipation); k_p→1 → all plastic.
	                          * Dissipation requires k_p > 1 - 2/(n+1), e.g. for n=1
	                          * any k_p > 0 gives dissipation. Default: 0.5 */

	/* ---- friction ---- */
	doublereal m_dMuS;

	/* ---- nodes ---- */
	unsigned                        m_nNodes;
	std::vector<const StructNode *> m_Nodes;

	/* ---- per-node plastic state ---- */
	std::vector<doublereal> m_dZmax;        /* committed max sinkage [m] */
	std::vector<doublereal> m_dZmax_trial;  /* trial (grows within timestep) */
	std::vector<doublereal> m_dKu;          /* unloading stiffness [N/m^(n+1)*m^n] */
	std::vector<doublereal> m_dZ0;          /* permanent sinkage offset [m] */

	/* ---- output ---- */
	mutable std::vector<Vec3> m_Forces;

	/* Bekker bearing coefficient */
	doublereal KBek(void) const { return m_dKc / m_dB + m_dKphi; }

	/* Recompute unloading parameters from a given z_max.
	 *
	 * Permanent sinkage:  z0 = k_plastic * z_max
	 * Unloading stiffness (ensures force continuity at z_max):
	 *   Ku = F_BEK(z_max) / (z_max - z0)
	 *      = A*KBek*z_max^n / (z_max*(1 - k_p))
	 *      = A*KBek*z_max^(n-1) / (1 - k_p)
	 *
	 * Energy dissipated per cycle (positive for k_p > 1 - 2/(n+1)):
	 *   E = A*KBek*z_max^(n+1) * [1/(n+1) - (1-k_p)/2]
	 */
	void UpdateUnloading(unsigned i, doublereal z_max) {
		doublereal z_ref = std::max(z_max, Z_FLOOR);
		m_dKu[i] = KBek() * std::pow(z_ref, m_dNbek - 1.0)
		         / (1.0 - m_dKplastic);
		m_dZ0[i] = m_dKplastic * z_max;
	}

public:
	SurfaceImpactElem(unsigned uLabel, const DofOwner *pDO,
		DataManager *pDM, MBDynParser& HP);
	virtual ~SurfaceImpactElem(void);

	virtual Elem::Type GetElemType(void) const { return Elem::LOADABLE; }
	virtual void Output(OutputHandler& OH) const;

	virtual void WorkSpaceDim(integer *piNumRows, integer *piNumCols) const {
		*piNumRows = 3 * m_nNodes;
		*piNumCols = 3 * m_nNodes;
	}

	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Commit plastic state after each converged timestep */
	virtual void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);

	virtual unsigned int iGetNumPrivData(void) const { return 0; }
	virtual int iGetNumConnectedNodes(void) const {
		return static_cast<int>(m_nNodes);
	}
	virtual void GetConnectedNodes(
		std::vector<const Node *>& connectedNodes) const
	{
		connectedNodes.resize(m_nNodes);
		for (unsigned i = 0; i < m_nNodes; i++) {
			connectedNodes[i] = m_Nodes[i];
		}
	}
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph) { }
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const { return 0; }
	virtual void InitialWorkSpaceDim(
		integer *piNumRows, integer *piNumCols) const
	{ *piNumRows = 0; *piNumCols = 0; }
	VariableSubMatrixHandler& InitialAssJac(
		VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
	{ WorkMat.SetNullMatrix(); return WorkMat; }
	SubVectorHandler& InitialAssRes(
		SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
	{ WorkVec.ResizeReset(0); return WorkVec; }
};

/* ========================================================================= */
/* Constructor                                                                 */
/* ========================================================================= */
SurfaceImpactElem::SurfaceImpactElem(unsigned uLabel, const DofOwner *pDO,
	DataManager *pDM, MBDynParser& HP)
: UserDefinedElem(uLabel, pDO),
m_Normal(Zero3), m_PlanePoint(Zero3),
m_dKc(0.), m_dKphi(0.), m_dNbek(1.), m_dB(1.), m_dANOR(0.),
m_dKplastic(0.5),
m_dMuS(0.), m_nNodes(0)
{
	if (HP.IsKeyWord("help")) {
		silent_cerr(
			"SurfaceImpactElem:\n"
			"  user defined: <label>, surface impact,\n"
			"      plane normal, <nx>, <ny>, <nz>,\n"
			"      plane point,  <px>, <py>, <pz>,\n"
			"      kc,         <kc>,\n"
			"      kphi,       <kphi>,\n"
			"      n,          <n>,\n"
			"      pad radius, <b>,\n"
			"      [ , area, <anor>, ]\n"
			"      [ , plastic ratio, <k_p>, ]\n"
			"      [ , friction, <mu_s>, ]\n"
			"      nodes, <N>, <label_1>, ..., <label_N>;\n"
			<< std::endl);
		if (!HP.IsArg()) throw NoErr(MBDYN_EXCEPT_ARGS);
	}

	/* plane normal */
	if (!HP.IsKeyWord("plane" "normal")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"plane normal\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_Normal = HP.GetVec3();
	doublereal dNorm = m_Normal.Norm();
	if (dNorm < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): zero-length normal at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_Normal *= (1.0 / dNorm);

	/* plane point */
	if (!HP.IsKeyWord("plane" "point")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"plane point\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_PlanePoint = HP.GetVec3();

	/* kc */
	if (!HP.IsKeyWord("kc")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"kc\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dKc = HP.GetReal();
	if (m_dKc < 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): kc must be >= 0 at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* kphi */
	if (!HP.IsKeyWord("kphi")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"kphi\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dKphi = HP.GetReal();
	if (m_dKphi <= 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): kphi must be > 0 at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* n */
	if (!HP.IsKeyWord("n")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"n\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dNbek = HP.GetReal();
	if (m_dNbek < 0.5) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): n must be >= 0.5 at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* pad radius */
	if (!HP.IsKeyWord("pad" "radius")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"pad radius\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dB = HP.GetReal();
	if (m_dB <= 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): pad radius must be > 0 at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* area (optional, default pi*b^2) */
	m_dANOR = M_PI_L * m_dB * m_dB;
	if (HP.IsKeyWord("area")) {
		m_dANOR = HP.GetReal();
		if (m_dANOR <= 0.) {
			silent_cerr("SurfaceImpactElem(" << uLabel
				<< "): area must be > 0 at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* plastic ratio k_p (optional, default 0.5)
	 * z0 = k_p * z_max  — fraction of max sinkage that remains permanently.
	 * Must satisfy k_p > 1 - 2/(n+1) for energy dissipation.
	 * For n=1: any k_p > 0 dissipates energy. */
	if (HP.IsKeyWord("plastic" "ratio")) {
		m_dKplastic = HP.GetReal();
		if (m_dKplastic <= 0. || m_dKplastic >= 1.) {
			silent_cerr("SurfaceImpactElem(" << uLabel
				<< "): plastic ratio must be in (0, 1) at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* friction (optional) */
	if (HP.IsKeyWord("friction")) {
		m_dMuS = HP.GetReal();
		if (m_dMuS < 0.) {
			silent_cerr("SurfaceImpactElem(" << uLabel
				<< "): friction must be >= 0 at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* node list */
	if (!HP.IsKeyWord("nodes")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"nodes\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	int nNodes = HP.GetInt();
	if (nNodes <= 0) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): number of nodes must be > 0 at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_nNodes = static_cast<unsigned>(nNodes);
	m_Nodes.resize(m_nNodes);
	m_Forces.resize(m_nNodes, Zero3);
	m_dZmax.resize(m_nNodes, 0.);
	m_dZmax_trial.resize(m_nNodes, 0.);
	m_dKu.resize(m_nNodes, 0.);
	m_dZ0.resize(m_nNodes, 0.);

	for (unsigned i = 0; i < m_nNodes; i++) {
		unsigned uNodeLabel = HP.GetInt();
		m_Nodes[i] = dynamic_cast<const StructNode *>(
			pDM->pFindNode(Node::STRUCTURAL, uNodeLabel));
		if (!m_Nodes[i]) {
			silent_cerr("SurfaceImpactElem(" << uLabel
				<< "): structural node " << uNodeLabel
				<< " not found at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		UpdateUnloading(i, 0.);
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

SurfaceImpactElem::~SurfaceImpactElem(void) { NO_OP; }

/* ========================================================================= */
/* AssRes                                                                      */
/* ========================================================================= */
SubVectorHandler& SurfaceImpactElem::AssRes(
	SubVectorHandler& WorkVec,
	doublereal /* dCoef */,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */)
{
	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iRow = 1;
	for (unsigned i = 0; i < m_nNodes; i++) {
		const StructNode *pNode = m_Nodes[i];
		integer iMomIdx = pNode->iGetFirstMomentumIndex();

		for (integer j = 1; j <= 3; j++) {
			WorkVec.PutRowIndex(iRow + j - 1, iMomIdx + j);
		}

		m_Forces[i] = Zero3;

		const Vec3& x = pNode->GetXCurr();
		const Vec3& v = pNode->GetVCurr();

		/* signed distance: positive = above plane (no contact) */
		doublereal d = (x - m_PlanePoint).Dot(m_Normal);
		if (d >= 0.) {
			iRow += 3;
			continue;
		}

		doublereal z = -d;   /* penetration depth >= 0 */

		/* ---- Bekker-Wong (mandatory, with plasticity) ---- */

		/* Update trial z_max monotonically within the timestep.
		 * Uses committed m_dZmax to decide the branch — trial only grows. */
		if (z > m_dZmax_trial[i]) {
			m_dZmax_trial[i] = z;
		}

		doublereal F_BEK, dFBEK_dz;

		if (z >= m_dZmax[i]) {
			/* Loading branch: virgin or deeper penetration */
			doublereal z_safe = std::max(z, Z_FLOOR);
			F_BEK    = m_dANOR * KBek() * std::pow(z_safe, m_dNbek);
			dFBEK_dz = m_dANOR * KBek() * m_dNbek
			         * std::pow(z_safe, m_dNbek - 1.0);
		} else {
			/* Unloading/reloading branch: elastic from committed z_max.
			 * Free zone z < z0: no force, no stiffness (zero Jacobian). */
			doublereal dz = z - m_dZ0[i];
			if (dz > 0.) {
				F_BEK    = m_dANOR * m_dKu[i] * dz;
				dFBEK_dz = m_dANOR * m_dKu[i];
			} else {
				F_BEK    = 0.;
				dFBEK_dz = 0.;   /* no contact in free zone */
			}
		}

		/* ---- Total normal force ---- */
		doublereal FN = F_BEK;

		(void)dFBEK_dz;   /* computed fresh in AssJac */

		Vec3 F_total = m_Normal * FN;

		/* ---- Coulomb friction on total FN ---- */
		if (m_dMuS > 0. && FN > 0.) {
			doublereal vdotn = v.Dot(m_Normal);
			Vec3 VT_vec = v - m_Normal * vdotn;
			doublereal VT_mag = VT_vec.Norm();
			if (VT_mag > 1.e-12) {
				F_total += VT_vec * (-m_dMuS * FN / VT_mag);
			}
		}

		m_Forces[i] = F_total;
		WorkVec.Add(iRow, F_total);
		iRow += 3;
	}

	return WorkVec;
}

/* ========================================================================= */
/* AssJac                                                                      */
/* ========================================================================= */
/*
 * Jacobian: J_jk = Jnn * n_j * n_k   (rows=momentum, cols=position)
 *   Loading:         Jnn = -dCoef * ANOR*KBek*n*z^(n-1)
 *   Unloading z>z0:  Jnn = -dCoef * ANOR*Ku
 *   Free zone z<z0:  Jnn = 0  (no contact force)
 * Friction Jacobian omitted.
 */
VariableSubMatrixHandler& SurfaceImpactElem::AssJac(
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */)
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	integer iNumRows, iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iRow = 1;
	for (unsigned i = 0; i < m_nNodes; i++) {
		const StructNode *pNode = m_Nodes[i];
		integer iMomIdx = pNode->iGetFirstMomentumIndex();
		integer iPosIdx = pNode->iGetFirstPositionIndex();

		for (integer j = 1; j <= 3; j++) {
			WM.PutRowIndex(iRow + j - 1, iMomIdx + j);
			WM.PutColIndex(iRow + j - 1, iPosIdx + j);
		}

		const Vec3& x = pNode->GetXCurr();

		doublereal d = (x - m_PlanePoint).Dot(m_Normal);
		if (d >= 0.) {
			iRow += 3;
			continue;
		}

		doublereal z = -d;

		/* Bekker stiffness (position-dependent, with dCoef).
		 * Loading:    dF/dz = ANOR*KBek*n*z^(n-1)
		 * Unloading:  dF/dz = ANOR*Ku   (only for z > z0; zero in free zone) */
		doublereal dFBEK_dz;
		if (z >= m_dZmax[i]) {
			doublereal z_safe = std::max(z, Z_FLOOR);
			dFBEK_dz = m_dANOR * KBek()
				* m_dNbek * std::pow(z_safe, m_dNbek - 1.0);
		} else if (z > m_dZ0[i]) {
			dFBEK_dz = m_dANOR * m_dKu[i];
		} else {
			dFBEK_dz = 0.;   /* free zone: no force gradient */
		}

		/* z = -d = -(x-P0)·n  →  dz/dx_k = -n_k
		 * dF_j/dx_k = dF/dz * (-n_k) * n_j  →  J_jk = -dCoef * dFBEK_dz * n_j * n_k */
		doublereal Jnn = -dCoef * dFBEK_dz;

		for (integer j = 1; j <= 3; j++) {
			for (integer k = 1; k <= 3; k++) {
				WM.IncCoef(iRow + j - 1, iRow + k - 1,
					Jnn * m_Normal(j) * m_Normal(k));
			}
		}

		iRow += 3;
	}

	return WorkMat;
}

/* ========================================================================= */
/* AfterConvergence — commit plastic state per node                           */
/* ========================================================================= */
void SurfaceImpactElem::AfterConvergence(
	const VectorHandler& /* X */,
	const VectorHandler& /* XP */)
{
	for (unsigned i = 0; i < m_nNodes; i++) {
		if (m_dZmax_trial[i] > m_dZmax[i]) {
			m_dZmax[i] = m_dZmax_trial[i];
			UpdateUnloading(i, m_dZmax[i]);
		}
		/* Reset trial to committed value for next timestep */
		m_dZmax_trial[i] = m_dZmax[i];
	}
}

/* ========================================================================= */
/* Output                                                                      */
/* ========================================================================= */
void SurfaceImpactElem::Output(OutputHandler& OH) const
{
	if (bToBeOutput() && OH.IsOpen(OutputHandler::LOADABLE)) {
		std::ostream& out = OH.Loadable();
		for (unsigned i = 0; i < m_nNodes; i++) {
			out << GetLabel()
				<< " " << m_Nodes[i]->GetLabel()
				<< " " << m_Forces[i]
				<< std::endl;
		}
	}
}

/* ========================================================================= */
/* Restart                                                                     */
/* ========================================================================= */
std::ostream& SurfaceImpactElem::Restart(std::ostream& out) const
{
	out << "user defined: " << GetLabel() << ", surface impact"
		<< ", plane normal, " << m_Normal
		<< ", plane point, "  << m_PlanePoint
		<< ", kc, "           << m_dKc
		<< ", kphi, "         << m_dKphi
		<< ", n, "            << m_dNbek
		<< ", pad radius, "   << m_dB
		<< ", area, "         << m_dANOR
		<< ", plastic ratio, "<< m_dKplastic
		<< ", friction, " << m_dMuS
		<< ", nodes, " << m_nNodes;
	for (unsigned i = 0; i < m_nNodes; i++) {
		out << ", " << m_Nodes[i]->GetLabel();
	}
	out << ";" << std::endl;
	return out;
}

/* ========================================================================= */
/* Reader and registration                                                     */
/* ========================================================================= */
struct SurfaceImpactElemRead : public UserDefinedElemRead {
	virtual UserDefinedElem *
	Read(unsigned uLabel, const DofOwner *pDO,
		DataManager *pDM, MBDynParser& HP) const
	{
		return new SurfaceImpactElem(uLabel, pDO, pDM, HP);
	}
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new SurfaceImpactElemRead;
	if (!SetUDE("surface" "impact", rf)) {
		delete rf;
		silent_cerr("SurfaceImpactElem: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}
	return 0;
}
