# How to contribute 
You can contribute to this project in many ways: 

* [contributing new features or improving existing ones;](## Development)  
* [subscribing and partecipating the pw-forum mailing list;](http://www.qe-forge.org/mailman/listinfo/pw_forum) 
* [reporting bug and proposing changes in the Issue section of the gitlab repository;](##Issues)
* [preparing new tests for the test suite](## Adding tests)

## Development 
If you want to contribute serious and non-trivial stuff ( or even simple and trivial  stuff ) you just have  to *fork* this repository; keep it updated; when your contribution is ready submit a merge request to this repository. At the moment the merge requests has to  be addressed to the master branch; this shall  change in the next future. 
     A basic guide on how to work with `git` can be found [here](https://docs.gitlab.com/ce/gitlab-basics/README.html). A  more thorough introduction to `git` is provided by [proGit](https://git-scm.com/book/en/v2) online e-book 

#### Proposed workflow 

   - register on [gitlab](https://gitlab.com/users/sign_in); 
   - fork the QEF/q-e project [?](https://docs.gitlab.com/ce/gitlab-basics/fork-project.html);
   - create a local copy of __your__ repository on your workstation(s). [?](https://docs.gitlab.com/ce/gitlab-basics/command-line-commands.html#clone-your-project);
   - create a branch where you work at your contribution[?](https://docs.gitlab.com/ce/gitlab-basics/start-using-git.html#create-a-branch); 
   - once your contribution is ready, if you haven't already done that, push it into __your__ repository [?](https://docs.gitlab.com/ce/gitlab-basics/start-using-git.html#send-changes-to-gitlab-com); 
   - submit a merge request to QEF/q-e[?](https://docs.gitlab.com/ce/gitlab-basics/add-merge-request.html#how-to-create-a-merge-request)
   
#### How to keep  updated  your repository

To be  sure that your work may  merged directlly to QEF\q-e repository it is important to keep it aligned to the changes occurring in the main repository. This is easier if you

* add QEF/q-e as a remote of your repository _e.g._ 

```
git remote add QEFqe https://gitlab.com/QEF/q-e.git 
git fetch upstream 
``` 
   
* create a branch which tracks QEF/q-e master *e.g.* 

```
git checkout -b QEFqe QEFqe/master
```

* before merging 
    - checkout to the tracking branch  
        ```
        git checkout QEFqe
        ```              
    - create a new branch identical to QEFqe and checkout to it   
        ```
        git checkout -b merge_code_contribution
        ``` 
    - merge your *code_contribution* branch in *merge_code_contribution*                  
        ```
        git merge code_contribution
        ```
    - if the merge has been successfull you can push your new branch to gitlab  
       ```
       git push origin merge_code_contribution
       ```      
    - you can now create your merge request from gitlab web interface 