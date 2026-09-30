#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include<streambuf>

/************************************************
 *  unit test of pseudo
 ***********************************************/

/**
 * - Tested Functions:
 *   - pseudo
 *   - ~pseudo
 *   - complete_default_h
 *   - complete_default_atom
 *   - complete_default_vl
 *   - complete_default
 *   - print_pseudo_h
 *   - print_pseudo_atom
 *   - print_pseudo_vl
 *   - print_pseudo
 *   - group_proj_by_l
 */

#define private public
#include "module_cell/read_pp.h"
#include "module_cell/atom_pseudo.h"
#undef private
class NCPPTest : public testing::Test
{
protected:
	std::unique_ptr<Pseudopot_upf> upf{new Pseudopot_upf};
	std::unique_ptr<Atom_pseudo> ncpp{new Atom_pseudo};
};

TEST_F(NCPPTest, SetPseudoH)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	PARAM.input.pseudo_rcut = 15.0;
	upf->read_pseudo_upf201(ifs, *ncpp);
	//set_pseudo_h
	upf->complete_default_h(*ncpp);

	if(!ncpp->has_so)
	{
		for (int i=0;i<ncpp->nchi;i++)
		{
			EXPECT_EQ(ncpp->nn[i],0);
			EXPECT_EQ(ncpp->jchi[i],0);
		}
		for (int i=0;i<ncpp->nbeta;i++)
		{
			EXPECT_EQ(ncpp->jjj[i],0);
		}
	}
	ifs.close();
}

TEST_F(NCPPTest, SetPseudoAtom)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	PARAM.input.pseudo_rcut = 15.0;
	upf->read_pseudo_upf201(ifs, *ncpp);
	//set_pseudo_atom
	upf->complete_default_h(*ncpp);
	upf->complete_default_atom(*ncpp);
	EXPECT_EQ(ncpp->rcut,PARAM.input.pseudo_rcut);

	if(!ncpp->nlcc)
	{
		for(int i=0;i<ncpp->mesh;i++)
		{
			EXPECT_EQ(ncpp->rho_atc[i],0.0);
		}
	}
	EXPECT_EQ(ncpp->msh,ncpp->mesh);
	ifs.close();
}

TEST_F(NCPPTest, SetPseudoNC)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	PARAM.input.pseudo_rcut = 15.0;
	// set pseudo nbeta = 0
	upf->read_pseudo_upf201(ifs, *ncpp);
	ncpp->nbeta = 0;
	upf->complete_default(*ncpp);
	EXPECT_EQ(ncpp->nh,0);
    // set pseudo nbeta > 0
	upf->read_pseudo_upf201(ifs, *ncpp);
    upf->complete_default(*ncpp);
	EXPECT_EQ(ncpp->nh,14);
	EXPECT_EQ(ncpp->kkbeta,132);
	ifs.close();
	
}

TEST_F(NCPPTest, PrintNC)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	PARAM.input.pseudo_rcut = 15.0;
	upf->read_pseudo_upf201(ifs, *ncpp);
    upf->complete_default(*ncpp);
    ifs.close();
	//print
	std::ofstream ofs;
	ofs.open("./tmp_log");
	ncpp->print_pseudo(ofs);
	ofs.close();
	ifs.open("./tmp_log");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str,testing::HasSubstr("psd  C"));
	EXPECT_THAT(str,testing::HasSubstr("pp_type  NC"));
	EXPECT_THAT(str,testing::HasSubstr("dft  PBE"));
	EXPECT_THAT(str,testing::HasSubstr("zv       4"));
	EXPECT_THAT(str,testing::HasSubstr("nchi 3"));
	EXPECT_THAT(str,testing::HasSubstr("nbeta    6"));
	EXPECT_THAT(str,testing::HasSubstr("dion :  nr=6 nc=6"));
	EXPECT_THAT(str,testing::HasSubstr("msh\t1247"));
	ifs.close();
	remove("./tmp_log");
}

TEST(NCPPStandaloneTest, GroupProjByL)
{
	// a pseudopotential that lists its projectors as l = 0,0,1,2,3,1, i.e. not grouped
	// by angular momentum (UPF does not require the grouping)
	pseudo pp;
	pp.nbeta = 6;
	pp.mesh = 4;
	pp.lll = {0, 0, 1, 2, 3, 1};
	pp.betar.create(6, 4);
	pp.dion.create(6, 6);
	for (int i = 0; i < 6; ++i)
	{
		for (int ir = 0; ir < 4; ++ir)
		{
			pp.betar(i, ir) = 10.0 * i + ir;
		}
		for (int j = 0; j < 6; ++j)
		{
			pp.dion(i, j) = 100.0 * i + j;
		}
	}

	pp.group_proj_by_l();

	// the projectors and the matching rows/columns of dion are permuted together
	const std::vector<int> expected_l = {0, 0, 1, 1, 2, 3};
	const std::vector<int> expected_source = {0, 1, 2, 5, 3, 4};
	for (int i = 0; i < 6; ++i)
	{
		EXPECT_EQ(pp.lll[i], expected_l[i]);
		for (int ir = 0; ir < 4; ++ir)
		{
			EXPECT_DOUBLE_EQ(pp.betar(i, ir), 10.0 * expected_source[i] + ir);
		}
		for (int j = 0; j < 6; ++j)
		{
			EXPECT_DOUBLE_EQ(pp.dion(i, j), 100.0 * expected_source[i] + expected_source[j]);
		}
	}

	// a pseudopotential that is already grouped is left untouched
	pp.group_proj_by_l();
	for (int i = 0; i < 6; ++i)
	{
		EXPECT_EQ(pp.lll[i], expected_l[i]);
		EXPECT_DOUBLE_EQ(pp.betar(i, 0), 10.0 * expected_source[i]);
	}
}
