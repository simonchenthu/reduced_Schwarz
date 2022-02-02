# A data-driven reduced order Schwarz method for nonlinear mutitscale PDEs

The source code for the paper [S. Chen, Q. Li, J. Lu, S. J. Wright. Manifold Learning and Nonlinear Homogenization. arXiv:2011.00568](https://arxiv.org/abs/2011.00568).

## Organization

- `src`: contains the local PDE solvers, the reduced order Schwarz solvers and all the sub routines
- `examples`: contains one example of semilinear elliptic equations and one example of nonlinear radiative transfer equations

The run time to generate offline dictionaries could be between several minutes to several hours depending on the parameters, e.g., the discretization and the dictionary size.

## Instructions for use

The instructions for running each case are as follows.

### Semilinear elliptic equations

1. Run elliptic_dictionary.m, which will generate dictionaries on each local patch and save them in data_elliptic
2. Run elliptic_ref.m, which will generate a reference solution and save it in data_elliptic
3. Run elliptic_online.m to solve semilinear elliptic equations by the reduced order Schwarz iteration

### Nonlinear radiative transfer equations

1. Run RTE_dictionary_interior.m and RTE_dictionary_boundary.m, which will generate dictionaries and save them in data_RTE
2. Run RTE_ref.m, which will generate a reference solution
3. Run RTE_online.m to solve the nonlinear radiative transfer equations by the reduced order Schwarz iteration

## Cite this work

If you use this code for academic research, you are encouraged to cite the following paper:

```
@article{ChLiLuWr:2020manifold,
  title={Manifold Learning and Nonlinear Homogenization},
  author={Chen, Shi and Li, Qin and Lu, Jianfeng and Wright, Stephen J},
  journal={arXiv preprint arXiv:2011.00568},
  year={2020}
}
```

## Questions

To get help on how to use the data or code, simply open an issue in the GitHub "Issues" section.