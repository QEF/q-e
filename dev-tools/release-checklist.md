## Release checklist for v.6.4

1. run all examples; ideally also check for discrepancies, or at least for crashes
2. verify that all README, README.md, etc. files contain updated information on the content of the relative package and directory
3. verify that all documentation files, in particular Doc/developer-man.tex and all user_guide.tex files, contain updated information. In particular, verify that there are no reference to removed or obsolete software and no missing references to new or changed software.
4. update the release number in developer_man.tex and in all user_guide.tex and other documentation that contains reference to version number
5. verify that input documentation (files INPUT_*.def) is updated
6. update Doc/release-notes with the release number and with updated information on what is new, changed, removed, etc.
7. Re-generate new documentation with "make doc"
8. verify that install/configure is updated and aligned with install/configure.ac
9. update version number in Modules/version.f90
10. set a git branch "qe-x.y[.z]" for version x.y[.z]
11. align master to develop, github to gitlab
12. make packages on gitlab and github
13. if there are changes to che schema, copy the new schema to 
quantumespresso@qe.safevps.it:/storage/vhosts/quantum-espresso.org/ns/qes
14. update the web site: add a piece of news, update pages Downloads, Roadmap, and any other page that needs to be updated, copy the new documentation to directory quantumespresso@qe.safevps.it:/storage/vhosts/quantum-espresso.org/htdocs/Doc
15. send a message to the mailing list, post to twitter, facebook, and whatnot
