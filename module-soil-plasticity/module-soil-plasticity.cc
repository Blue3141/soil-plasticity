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
 * Surface impact element for MBDyn.
 *
 * Models the normal reaction force of a planar surface on a set of structural
 * nodes using a hydrodynamic (quadratic velocity) formulation:
 *
 *   FNH = ACCA2 * ANOR * VN^2 * tanh(VN / VN0)
 *
 * where:
 *   FNH   = normal reaction force [N]
 *   ACCA2 = 0.5 * rho * Cd  (rho = fluid/soil density, Cd = drag coefficient;
 *            for a flat plate orthogonal to flow: Cd = 1.2)
 *   ANOR  = effective contact area [m^2]
 *   VN    = approach velocity normal to the plane (positive = into surface) [m/s]
 *   VN0   = reference velocity for tanh stabilisation (typically 10% of max
 *            expected VN) [m/s]
 *
 * The tanh factor acts as a smooth sign function: for VN >> VN0 the formula
 * reduces to the standard quadratic drag law FNH = ACCA2*ANOR*VN*|VN|.
 * For VN < 0 (node moving away from surface) tanh is negative and the force
 * tends to zero — no adhesion.
 *
 * Optional static friction resists tangential slip:
 *   F_fric = -mu_s * FNH * (VT / |VT|)
 *
 * The element activates only when a node has penetrated the plane
 * (signed distance d = (x - P0).n < 0).
 *
 * Input syntax:
 *   user defined: <label>, surface impact,
 *       plane normal, <nx>, <ny>, <nz>,   # outward unit normal
 *       plane point,  <px>, <py>, <pz>,   # any point on the plane
 *       ACCA2,  <acca2>,                  # 0.5*rho*Cd [kg/m^3]
 *       area,   <anor>,                   # contact area per node [m^2]
 *       VN0,    <vn0>,                    # reference velocity [m/s]
 *       [ , friction, <mu_s>, ]           # static friction coeff [-]
 *       nodes, <N>,                       # number of nodes
 *           <label_1>, ..., <label_N>     # structural node labels
 *       ;
 *
 * In control data:
 *   loadable elements: 1;
 *
 * Registration key: "surface impact"
 */

#include "mbconfig.h"

#include <cmath>
#include <limits>
#include <vector>

#include "dataman.h"
#include "userelem.h"
#include "strnode.h"

/* ========================================================================= */
class SurfaceImpactElem : public UserDefinedElem {
private:
	Vec3       m_Normal;      /* outward unit normal of the plane              */
	Vec3       m_PlanePoint;  /* a point on the plane                          */
	doublereal m_dACCA2;      /* 0.5 * rho * Cd [kg/m^3]                      */
	doublereal m_dANOR;       /* effective contact area [m^2]                  */
	doublereal m_dVN0;        /* reference velocity for tanh [m/s]             */
	doublereal m_dMuS;        /* static friction coefficient [-], 0 = no fric  */

	unsigned                         m_nNodes;
	std::vector<const StructNode *>  m_Nodes;

	/* Force on node i, computed in AssRes — stored for Output */
	mutable std::vector<Vec3> m_Forces;

public:
	SurfaceImpactElem(unsigned uLabel, const DofOwner *pDO,
		DataManager *pDM, MBDynParser& HP);
	virtual ~SurfaceImpactElem(void);

	/* ---- element type ---- */
	virtual Elem::Type GetElemType(void) const { return Elem::LOADABLE; }

	/* ---- output ---- */
	virtual void Output(OutputHandler& OH) const;

	/* ---- workspace sizing ---- */
	virtual void WorkSpaceDim(integer *piNumRows, integer *piNumCols) const {
		/* 3 translational DOFs per node, no couples */
		*piNumRows = 3 * m_nNodes;
		*piNumCols = 3 * m_nNodes;
	}

	/* ---- Jacobian ---- */
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* ---- residual ---- */
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* ---- private data (none) ---- */
	virtual unsigned int iGetNumPrivData(void) const { return 0; }

	/* ---- connected nodes ---- */
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

	/* ---- misc required virtuals ---- */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph) { }

	virtual std::ostream& Restart(std::ostream& out) const;

	/* ---- initial assembly (not used) ---- */
	virtual unsigned int iGetInitialNumDof(void) const { return 0; }
	virtual void InitialWorkSpaceDim(
		integer *piNumRows, integer *piNumCols) const
	{
		*piNumRows = 0;
		*piNumCols = 0;
	}
	VariableSubMatrixHandler& InitialAssJac(
		VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}
	SubVectorHandler& InitialAssRes(
		SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
	{
		WorkVec.ResizeReset(0);
		return WorkVec;
	}
};

/* ========================================================================= */
/* Constructor / parser                                                        */
/* ========================================================================= */
SurfaceImpactElem::SurfaceImpactElem(unsigned uLabel, const DofOwner *pDO,
	DataManager *pDM, MBDynParser& HP)
: UserDefinedElem(uLabel, pDO),
m_Normal(Zero3), m_PlanePoint(Zero3),
m_dACCA2(0.), m_dANOR(0.), m_dVN0(1.), m_dMuS(0.),
m_nNodes(0)
{
	if (HP.IsKeyWord("help")) {
		silent_cerr(
			"SurfaceImpactElem:\n"
			"  user defined: <label>, surface impact,\n"
			"      plane normal, <nx>, <ny>, <nz>,\n"
			"      plane point,  <px>, <py>, <pz>,\n"
			"      ACCA2,  <acca2>,     # 0.5*rho*Cd [kg/m^3]\n"
			"      area,   <anor>,      # contact area [m^2]\n"
			"      VN0,    <vn0>,       # reference velocity [m/s]\n"
			"      [ , friction, <mu_s>, ]\n"
			"      nodes, <N>,\n"
			"          <label_1>, ..., <label_N>;\n"
			<< std::endl);
		if (!HP.IsArg()) {
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	/* -- plane normal -- */
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
			<< "): zero-length normal vector at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_Normal *= (1.0 / dNorm);   /* normalise */

	/* -- plane point -- */
	if (!HP.IsKeyWord("plane" "point")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"plane point\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_PlanePoint = HP.GetVec3();

	/* -- ACCA2 -- */
	if (!HP.IsKeyWord("ACCA2")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"ACCA2\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dACCA2 = HP.GetReal();
	if (m_dACCA2 <= 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): ACCA2 must be positive at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* -- area -- */
	if (!HP.IsKeyWord("area")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"area\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dANOR = HP.GetReal();
	if (m_dANOR <= 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): area must be positive at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* -- VN0 -- */
	if (!HP.IsKeyWord("VN0")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"VN0\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_dVN0 = HP.GetReal();
	if (m_dVN0 <= 0.) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): VN0 must be positive at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* -- friction (optional) -- */
	m_dMuS = 0.;
	if (HP.IsKeyWord("friction")) {
		m_dMuS = HP.GetReal();
		if (m_dMuS < 0.) {
			silent_cerr("SurfaceImpactElem(" << uLabel
				<< "): friction must be >= 0 at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* -- node list -- */
	if (!HP.IsKeyWord("nodes")) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): \"nodes\" expected at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	int nNodes = HP.GetInt();
	if (nNodes <= 0) {
		silent_cerr("SurfaceImpactElem(" << uLabel
			<< "): number of nodes must be positive at line "
			<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_nNodes = static_cast<unsigned>(nNodes);
	m_Nodes.resize(m_nNodes);
	m_Forces.resize(m_nNodes, Zero3);

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
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

SurfaceImpactElem::~SurfaceImpactElem(void)
{
	NO_OP;
}

/* ========================================================================= */
/* AssRes — residual contribution                                              */
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

		/* row indices for the 3 translational momentum equations */
		for (integer j = 1; j <= 3; j++) {
			WorkVec.PutRowIndex(iRow + j - 1, iMomIdx + j);
		}

		m_Forces[i] = Zero3;   /* default: no contact */

		const Vec3& x = pNode->GetXCurr();
		const Vec3& v = pNode->GetVCurr();

		/* signed distance from plane: positive = above (no contact) */
		doublereal d = (x - m_PlanePoint).Dot(m_Normal);
		if (d >= 0.) {
			iRow += 3;
			continue;
		}

		/* approach velocity: positive = moving into the surface */
		doublereal VN = -(v.Dot(m_Normal));

		/* normal force magnitude — tanh smooths VN*|VN| around zero */
		doublereal th  = std::tanh(VN / m_dVN0);
		doublereal FNH = m_dACCA2 * m_dANOR * VN * VN * th;

		/* no adhesion: clamp negative force to zero */
		if (FNH < 0.) {
			FNH = 0.;
		}

		/* normal force vector: pushes node away from surface (+n direction) */
		Vec3 F_total = m_Normal * FNH;

		/* tangential friction */
		if (m_dMuS > 0. && FNH > 0.) {
			doublereal vdotn = v.Dot(m_Normal);
			Vec3 VT_vec = v - m_Normal * vdotn;
			doublereal VT_mag = VT_vec.Norm();

			if (VT_mag > 1.e-12) {
				/* static friction opposes tangential motion */
				F_total += VT_vec * (-m_dMuS * FNH / VT_mag);
			}
		}

		m_Forces[i] = F_total;
		WorkVec.Add(iRow, F_total);
		iRow += 3;
	}

	return WorkVec;
}

/* ========================================================================= */
/* AssJac — Jacobian contribution                                              */
/* ========================================================================= */
/*
 * The normal force depends on the approach velocity VN = -(v·n).
 *
 * dFNH/dVN = ACCA2 * ANOR * VN * (2*tanh + VN/VN0 * (1 - tanh^2))
 *
 * Jacobian of the force vector w.r.t. velocity (3x3, per node):
 *   dF_normal/dv = -(dFNH/dVN) * (n ⊗ n)
 *
 * This is negative-semidefinite (damping-like), which is stabilising.
 * The friction Jacobian is omitted (approximation acceptable for convergence).
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
		const Vec3& v = pNode->GetVCurr();

		doublereal d = (x - m_PlanePoint).Dot(m_Normal);
		if (d >= 0.) {
			iRow += 3;
			continue;
		}

		doublereal VN    = -(v.Dot(m_Normal));
		doublereal th    = std::tanh(VN / m_dVN0);
		doublereal sech2 = 1.0 - th * th;

		/* dFNH/dVN: derivative of normal force magnitude w.r.t. approach velocity */
		doublereal dFdVN = m_dACCA2 * m_dANOR
			* VN * (2.0 * th + (VN / m_dVN0) * sech2);
		if (dFdVN < 0.) dFdVN = 0.;

		/*
		 * Velocity Jacobian (no dCoef factor):
		 *   dF_normal_j / dv_k = -(dFNH/dVN) * n_j * n_k
		 *
		 * Converted to position Jacobian via xp = (x - x_old)/dCoef:
		 *   dF/dx_k = (dF/dv_k) / dCoef
		 * In AssJac we add dCoef * dF/dx = dF/dv — so NO dCoef factor here.
		 */
		for (integer j = 1; j <= 3; j++) {
			for (integer k = 1; k <= 3; k++) {
				WM.IncCoef(iRow + j - 1, iRow + k - 1,
					-dFdVN * m_Normal(j) * m_Normal(k));
			}
		}

		iRow += 3;
	}

	return WorkMat;
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
		<< ", ACCA2, "  << m_dACCA2
		<< ", area, "   << m_dANOR
		<< ", VN0, "    << m_dVN0
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
