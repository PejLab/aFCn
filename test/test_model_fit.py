import unittest
import numpy as np
from afcn import model

# TODO: I haven't tested the ability to infer the correct
# parameters from data

class TestLinearized(unittest.TestCase):
    j_vars = 5
    n_samples = 1000
    q_values = j_vars * n_samples

    def setUp(self):
        if not hasattr(self, "rng"):
            self.rng = np.random.default_rng()

        self.hap_one = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.hap_two = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.effect_sizes = self.rng.normal(0,2,size=self.j_vars)

        self.y = (np.exp(self.hap_one @ self.effect_sizes)
                  + np.exp(self.hap_one @ self.effect_sizes))


    def test_output_types_properties(self):

        lsq = model._linear_expansion_model(self.hap_one,
                                            self.hap_two,
                                            self.y)

        self.assertIsInstance(lsq, dict)
        self.assertEqual(len(lsq), 2)
        for key in ["pars", "rank"]:
            self.assertTrue(key in lsq)


        self.assertEqual(lsq["pars"].size, self.j_vars+1)
        self.assertEqual(lsq["rank"], self.j_vars+1)

    # TODO
    def test_out_parameter_value(self):
        pass


class TestObj(unittest.TestCase):
    j_vars = 5
    n_samples = 1000
    q_values = j_vars * n_samples

    def setUp(self):
        if not hasattr(self, "rng"):
            self.rng = np.random.default_rng()

        self.hap_one = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.hap_two = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.effect_sizes = self.rng.normal(0,2,size=self.j_vars)

        self.y = (np.exp2(self.hap_one @ self.effect_sizes)
                  + np.exp2(self.hap_one @ self.effect_sizes))

    def test_variables_and_output_type(self):

        reg = "l1"
        regconst = 2

        input_vals = (self.hap_one, self.hap_two, self.y,
                      reg, 2)
        f = model._obj(*input_vals)
        self.assertTrue(callable(f))

        # test to make sure that the variables associated with
        # the closure are those that I intended.  Moreover,
        # the first elment of the closure is the lambda function
        # as reg and reg_const are not used in _g, they are not
        # stored in the closure.  Consequently, I am only testing
        # the matching of haplotypes and y

        for v, c in zip(input_vals, f.__closure__[1:]):

            self.assertTrue((v == c.cell_contents).all())
        

    def test_no_penalty(self):
        # Here we test whether a penalty term returns the intended
        # value.  This is done by picking output values that match
        # model predictions, i.e. the sum squared residuals (ssr) is
        # zero.  The value that is returned by the objective function
        # is then penalty term.
        reg = None
        regconst = None
        ssr_no_penalty = 0

        # recall the number of parameters is j_vars +1, where
        # the 1 accounts for the log reference abundance
        pars = np.ones(self.j_vars+1)

        # when all the haplotypes match the reference
        # allele, then I expect the model prediction to produce
        # 2* exp2(log_ref_const)  = 2 * exp2(1) = 4

        y = np.zeros(self.n_samples) + np.log2(1. + 4)

        haplotypes = np.zeros(shape=self.hap_one.shape)

        f = model._obj(haplotypes, haplotypes, y, reg, regconst)

        # Under no penalty, the returned value from the objective
        # should be ssr_no_penalty
        self.assertAlmostEqual(f(pars), ssr_no_penalty)
        

        # 2* 2^(log_ref_const)  = 2 * 2^(-1) = 1
        pars = -np.ones(self.j_vars + 1)
        y = np.zeros(self.n_samples) + np.log2(1 + 1)

        f = model._obj(haplotypes, haplotypes, y, reg, regconst)
        self.assertEqual(f(pars), ssr_no_penalty)


    def test_l1_penalty(self):
        # Here we test whether a penalty term returns the intended
        # value.  This is done by picking output values that match
        # model predictions, i.e. the sum squared residuals (ssr) is
        # zero.  The value that is returned by the objective function
        # is then penalty term.
        reg = "l1"
        regconst = 3.45

        # recall the number of parameters is j_vars +1, where
        # the 1 accounts for the log reference abundance
        pars = np.ones(self.j_vars+1)

        # when all the haplotypes are match the reference
        # allele, then I expect the model prediction to produce
        # 2* exp2(log_ref_const)  = 2 * exp2(1) = 4
        y = np.zeros(self.n_samples) + np.log2(1 + 4)

        haplotypes = np.zeros(shape=self.hap_one.shape)

        # given these inputs I expect that the only term 
        # of the objective that is non-zero
        # is the penalty term.  As each parameter = 1, then the 
        # penalty term under l1 is simply regconst * (j_vars + 1)
        f = model._obj(haplotypes, haplotypes, y, reg, regconst)

        self.assertEqual(f.__closure__[0].cell_contents.__closure__[0].cell_contents,
                         regconst)

        self.assertEqual(f(pars), regconst*(self.j_vars + 1))
        

        # 2* exp(log_ref_const)  = 2 * exp2(-1) = 1
        pars = -np.ones(self.j_vars+1)
        y = np.zeros(self.n_samples) + np.log2(1 + 1)

        f = model._obj(haplotypes, haplotypes, y, reg, regconst)
        self.assertEqual(f(pars), regconst*(self.j_vars + 1))


    def test_l2_penalty(self):
        # Here we test whether a penalty term returns the intended
        # value.  This is done by picking output values that match
        # model predictions, i.e. the sum squared residuals (ssr) is
        # zero.  The value that is returned by the objective function
        # is then penalty term.
        reg = "l2"
        regconst = 2.87

        # recall the number of parameters is j_vars +1, where
        # the 1 accounts for the log reference abundance
        pars = np.ones(self.j_vars+1)

        # when all the haplotypes are match the reference
        # allele, then I expect the model prediction to produce
        # 2* exp(log_ref_const)  = 2 * exp2(1) = 4
        y = np.zeros(self.n_samples) + np.log2(1 + 4)

        haplotypes = np.zeros(shape=self.hap_one.shape)

        # given these inputs I expect that the only nonzero
        # term of the objective is the penalty term
        # As each parameter = 1, then the 
        # penalty term under l2 is simply regconst * (j_vars + 1)
        f = model._obj(haplotypes, haplotypes, y, reg, regconst)

        self.assertEqual(f.__closure__[0].cell_contents.__closure__[0].cell_contents,
                         regconst)

        self.assertAlmostEqual(f(pars), regconst*(self.j_vars + 1))
        

        # 2* exp(log_ref_const)  = 2 * exp2(-1) = 1
        pars = -np.ones(self.j_vars+1)
        y = np.zeros(self.n_samples) + np.log2(1 + 1)

        f = model._obj(haplotypes, haplotypes, y, reg, regconst)
        self.assertAlmostEqual(f(pars), regconst*(self.j_vars + 1))


    def test_penalty_exception(self):
        # Incoherent penalty parameters should result in an exception being
        # raised
        reg = "l2"
        regconst = None

        with self.assertRaises(ValueError):
            f = model._obj(self.hap_one, self.hap_two, self.y, reg, regconst)

        with self.assertRaises(ValueError):
            f = model._obj(self.hap_one, self.hap_two, self.y, reg, "a")

        with self.assertRaises(ValueError):
            f = model._obj(self.hap_one, self.hap_two, self.y, reg, [1,2])

        with self.assertRaises(ValueError):
            f = model._obj(self.hap_one, self.hap_two, self.y, None, 1)


class TestFit(unittest.TestCase):
    j_vars = 5
    n_samples = 1000
    q_values = j_vars * n_samples

    def setUp(self):
        if not hasattr(self, "rng"):
            self.rng = np.random.default_rng()

        self.hap_one = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.hap_two = self.rng.choice([0,1],
                replace=True,
                size=self.q_values).reshape(self.n_samples, self.j_vars)

        self.effect_sizes = self.rng.normal(0,2,size=self.j_vars)

        self.y = (np.exp(self.hap_one @ self.effect_sizes)
                  + np.exp(self.hap_one @ self.effect_sizes))

    def test_fit_input_validation(self):
        with self.assertRaises(ValueError):
            self.hap_one[0,0] = 2

            model.fit(self.hap_one, self.hap_two, self.y) 

        with self.assertRaises(ValueError):
            self.hap_one[0,0] = -1

            model.fit(self.hap_one, self.hap_two, self.y) 

        with self.assertRaises(ValueError):
            self.hap_one = self.hap_one[:, :-1]

            model.fit(self.hap_one, self.hap_two, self.y) 


        with self.assertRaises(ValueError):
            self.y = self.y[:-1]

            model.fit(self.hap_one, self.hap_two, self.y) 


if __name__ == "__main__":
    unittest.main()
