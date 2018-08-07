# How to contribute
You can contribute to this project in many ways:

* [contributing new features or improving existing ones;](#development)
* [subscribing and partecipating to the users mailing list;](https://lists.quantum-espresso.org/mailman/listinfo/users)
* [reporting bugs and proposing changes in the Issue section of the gitlab repository;](#creating-issues)
* [preparing new tests for the test suite](#adding tests)

## Development
If you want to contribute serious and non-trivial stuff ( or even simple and trivial  stuff ) you just have  to *fork* this repository; keep it updated;
when your contribution is ready,  submit a merge request to the development branch  of this repository.
After some basic tests ran by gitlab CI and approval your changes will be merged in develop.
When  the whole  test-suite has been tested your contribution  will be merged  to the master branch.

A basic guide on how to work with `git` can be found [here](https://docs.gitlab.com/ce/gitlab-basics/README.html). A  more thorough introduction to `git` is provided by [proGit](https://git-scm.com/book/en/v2) online e-book



#### Proposed workflow

   - register on [gitlab](https://gitlab.com/users/sign_in);
   - [fork the QEF/q-e project](https://docs.gitlab.com/ce/gitlab-basics/fork-project.html);
   - [create a local copy of __your__ repository on your workstation(s)](https://docs.gitlab.com/ce/gitlab-basics/command-line-commands.html#clone-your-project);
   - [once your contribution is ready push it into __your__ repository](https://docs.gitlab.com/ce/gitlab-basics/start-using-git.html#send-changes-to-gitlab-com);
   - [create a  merge request to the `develop` branch of QEF/q-e](https://docs.gitlab.com/ce/gitlab-basics/add-merge-request.html#how-to-create-a-merge-request)


## Development tools
[Here](dev-tools/) you can find several tools that will assist you while contributing to the QE source code.

## Creating  Issues

You can report bugs and propose new developments posting on the [issue]( https://gitlab.com/QEF/q-e/issues)
section of this repository. Partecipation to the issue discussions are also a welcome contribution.
When reporting bugs try to help other partecipant to reproduce the problem by providing input and output files.
