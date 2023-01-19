# How to Contribute to ThermoEngine

We welcome contributions to this project and ask that you follow these procedures and [code of conduct](#code-of-conduct-for-contributors) when making your contribution.

### 1️⃣ Fork the repository.

The repository is publicly assessible, but we discourage branching or modfying the master branch, so the first step in contributing is to fork the repository to your private workspace on GitLab.
  1. Click the **Fork** button on the top right.
  1. Choose a namespace (your home on GitLab) from the list provided. This action creates a copy of the repository in that location.
  1. (Recommended) Rename your forked repository to, for example, *ThermoEngineMyInitials*, to avoid confusing it with the source repository.
     1. On GitLab, go to **Settings > General**, enter a new project name, and click **Save Changes**.
     1. (IMPORTANT) Make sure to keep the path name the same as your project name. Go to **Settings > General**, and expand the **Advanced** section. Under **Change path**, enter the path name, and click **Change path**.

### 2️⃣ Set up your forked repository to mirror the source repository.

**Note: The automatic mirroring setup described in this section works only if you have a paid GitLab membership. You can still work quite effectively without this step, so do not worry if you do not have a GitLab membership.**

After forking, set up your forked repository to mirror the source repository. This ensures that all future updates to the source repository are included in your fork. Avoid working directly in the master or Documentation branch of your forked repository to ensure that the mirroring process does not break or encounter a conflict as you develop code.

To set up mirroring on GitLab:
   1. At your forked repository, go to **Settings > Repository**, and expand the **Mirroring repositories** section.
   1. Enter https://gitlab.com/ENKI-portal/ThermoEngine in the **Git repository URL** field.
   1. From the **Mirror direction** dropdown menu, choose **Pull**.
   1. Click **Mirror repository**.

### 3️⃣ Restrict access if desired.

If you would like to restrict access to your project, it is important that you *do not reset access to Private*. That action will break the fork relationship. Instead, do the following procedure.

To restrict access to your project:
   1. Go to **Settings > General**, and expand the **Visibility, project, features, permissions** section.
   1. (IMPORTANT!) Leave **Project visibility** set to **Public**.
   1. Set permissions for **Issues, Repository, Forks, Pipelines, Wiki, Snippets and Pages** to **Only Project Members**.
   1. Click **Save changes**.

Your project files will be visible only to you.

### 4️⃣ Create a development branch from master in your forked repository.

Now that your mirror is set up, create a development branch from master in your forked repository, and give it a name, which we will assume is *develop*. Consider adopting a development workflow such as that recommended by [GitFlow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

### 5️⃣ Checkout the *develop* branch.

Check out the *develop* branch, and create your new code feature in that branch. Most users will want to work locally, rather than at GitLab.

To work locally:

   1. Clone the forked repository on GitLab to your local system using Git commands in a terminal or some GUI Git user interface (e.g., Tower, GitKrachen).
   1. Navigate to the local cloned source directory.
   1. Check out the *develop* branch of the clone.
   1. Revise code in the locally cloned workspace.

To setup a convenient local development environment, it is recommended that you setup automatic login for your local machine (using either a personal access token over https or ssh-keys over ssh). Both methods work, but we recommend the standard approach using personal access tokens.

To setup a automatic login using a personal access token:
   1. Login to GitLab.
   1. Click on your profile picture in the upper right and select *Edit profile*.
   1. Click on the *Access Tokens* button in the vertical left-hand menu bar.
   1. Select a Token name (e.g. 'work-laptop') and select the desired scope (*api* gives complete control for all possible use-cases).
   1. Click *Create personal access token*
   1. Copy the token value to a secure location (e.g. password manager). This token must be kept private as it grants full access to your repo without any passwords needed. *Do not lose the token as it is unrecoverable from GitLab. If you do lose it, you must use these steps to generate another.*
   1. Store your token in your local git repo using the following command at the command line:
      `git remote set-url origin https://oauth2:TOKEN@gitlab.com/USERNAME/YOUR_REPO.git`
      - be sure to replace TOKEN, USERNAME, and YOUR_REPO with your personal info.


### 6️⃣ After finishing your work, add it back to the repository.

When you are ready, add your code additions or revisions back to the ThermoEngine repository.

To add your code additions or revisions back to the ThermoEngine repository:
1. Do one of the following:
   - If working locally:  
     1. Use the Git **pull** command to update the cloned master branch from the forked repository master branch.
     1. Use Git to **merge** your cloned master branch into your cloned *develop* branch.  These steps ensure that the latest changes to master in the source repository are reflected in your code revisions.
     1. Use Git to **commit** the *develop* branch.
     1. Use the Git **push** command to send code revisions on that branch back up to the forked repository.
   - If working on GitLab:
     - Merge your forked repository master branch into your *develop* branch.  This step ensures that the latest changes to master in the source repository are reflected in your code revisions.   
1. In GitLab, navigate to the top view of your forked repository.
1. In the left-hand pane, click **Merge requests** and then click **New merge request**.
1. For **Source branch**, choose your forked project’s *develop* branch. For **Target branch**, choose the original Enki-portal/ThermoEngine master branch.
1. Click **Compare branches and continue**. You can now optionally add labels, add a milestone, and assign the merge request to someone who can review your changes.
1. Click **Submit merge request**.

### 7️⃣ Resolve any conflicts as instructed by maintainer.
At this stage, your merge request will be reviewed by a maintainer of the ThermoEngine repository, and one of the following will occur:
* The changes will be accepted.
* The maintainer will either resolve any conflicts or redirect the request back to you with instructions for resolution.

### 8️⃣ Congratulate yourself!

That's it! Congratulate yourself.  You have just made a contribution to the ThermoEngine code base and the ENKI project.


## Code of Conduct for Contributors

Please note that we have a code of conduct. Please follow it in all your interactions with the project.

### Our pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

### Our standards

Examples of behavior that contributes to creating a positive environment:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants:

* Using sexualized language or imagery, and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct that could reasonably be considered inappropriate in a professional setting

### Our responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project email address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by [emailing the project team](mailto:incoming+enki-portal-thermoengine-2313798-issue-@incoming.gitlab.com). All complaints will be reviewed and investigated and will result in a response deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version].

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
