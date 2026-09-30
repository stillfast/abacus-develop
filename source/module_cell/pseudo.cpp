#include "pseudo.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

#include <algorithm>
#include <vector>

pseudo::pseudo()
{
}

pseudo::~pseudo()
{
}

void pseudo::group_proj_by_l()
{
    if (nbeta <= 1 || static_cast<int>(lll.size()) < nbeta || betar.nr < nbeta || betar.nc < mesh || dion.nr < nbeta
        || dion.nc < nbeta)
    {
        return;
    }

    // order[i] is the projector that has to take the i-th place once the projectors are
    // grouped by l. std::stable_sort keeps the relative order of the projectors sharing
    // an l, which is the order of their zeta index.
    std::vector<int> order(nbeta);
    for (int i = 0; i < nbeta; ++i)
    {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&](const int i, const int j) { return lll[i] < lll[j]; });

    bool grouped = true;
    for (int i = 0; i < nbeta; ++i)
    {
        if (order[i] != i)
        {
            grouped = false;
            break;
        }
    }
    if (grouped)
    {
        return;
    }

    // The projectors of a fully relativistic pseudopotential are pairs (l,j), and the
    // augmentation charges of an ultrasoft one are stored per pair of projectors; both
    // would need to be permuted as well.
    if (has_so || tvanp)
    {
        ModuleBase::WARNING_QUIT("pseudo::group_proj_by_l",
                                 "the projectors of this pseudopotential are not grouped by angular momentum, which "
                                 "ABACUS only supports for norm-conserving pseudopotentials without spin-orbit "
                                 "coupling; please list the PP_BETA blocks by increasing angular momentum and "
                                 "permute PP_DIJ accordingly (see the '--sort-by-l' option of the SIAB tools)");
    }

    ModuleBase::matrix betar_new(nbeta, mesh);
    ModuleBase::matrix dion_new(nbeta, nbeta);
    std::vector<int> lll_new(nbeta);
    std::string lll_read;
    for (int i = 0; i < nbeta; ++i)
    {
        lll_new[i] = lll[order[i]];
        lll_read += std::to_string(lll[i]);
        if (i != nbeta - 1)
        {
            lll_read += ",";
        }
        for (int ir = 0; ir < mesh; ++ir)
        {
            betar_new(i, ir) = betar(order[i], ir);
        }
        for (int j = 0; j < nbeta; ++j)
        {
            dion_new(i, j) = dion(order[i], order[j]);
        }
    }
    betar = betar_new;
    dion = dion_new;
    lll = lll_new;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "note: projectors regrouped by angular momentum (read as l)",
                                lll_read);
}

void pseudo::print_pseudo(std::ofstream& ofs)
{
	print_pseudo_vl(ofs);
	ofs << "\n pseudo : ";
	ofs << "\n kkbeta	" << kkbeta;
	ofs << "\n nh  " << nh;
	output::printr1_d(ofs, " lll : ", lll.data(), nbeta);
	output::printrm(ofs, " betar : ", betar);
	output::printrm(ofs, " dion : ", dion);
	ofs << "\n ----------------------";
}

void pseudo::print_pseudo_atom(std::ofstream &ofs)
{
	print_pseudo_h(ofs);
	ofs << "\n pseudo_atom : ";
	ofs << "\n msh	" << msh;
//	ofs	<< "\n nchi	" << nchi;
	output::printr1_d(ofs, " r : ", r.data(), mesh);
	output::printr1_d(ofs, " rab : ", rab.data(), mesh);
	output::printr1_d(ofs, " rho_atc : ", rho_atc.data(), mesh);
	output::printr1_d(ofs, " rho_at : ", rho_at.data(), mesh);
	output::printr1_d(ofs," jchi : ", jchi.data(), nchi);
	output::printrm(ofs, " chi : ", chi);
	ofs << "\n ----------------------";
}


void pseudo::print_pseudo_vl(std::ofstream &ofs)
{
	ofs << "\n pseudo_vl:";
	print_pseudo_atom(ofs);
	output::printr1_d(ofs, "vloc_at : ", vloc_at.data(), mesh);
	ofs << "\n ----------------------------------- ";
}

void pseudo::print_pseudo_h(std::ofstream &ofs)
{
    ofs << "\n pseudo_info :";
    ofs << "\n nv       " << nv;
    ofs << "\n psd  " << psd;
    ofs << "\n pp_type  " << pp_type;
    ofs << "\n tvanp    " << tvanp;
    ofs << "\n nlcc " << nlcc;
    ofs << "\n dft  " << xc_func;
    ofs << "\n zv       " << zv;
    ofs << "\n etotps   " << etotps;
    ofs << "\n ecutwfc " << ecutwfc;
    ofs << "\n ecutrho " << ecutrho;
    ofs << "\n lmax " << lmax;
    ofs << "\n mesh " << mesh;
    ofs << "\n nchi " << nchi;
    ofs << "\n nbeta    " << nbeta;
//  out.printr1_d(ofs," els: ", els, nchi);
    output::printr1_d(ofs, " lchi: ", lchi.data(), nchi);
    output::printr1_d(ofs, " oc: ", oc.data(), nchi);
    ofs << "\n ----------------------";
}

