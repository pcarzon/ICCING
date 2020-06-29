This is a Tutorial for the use of Git Hub as retaining to the ICCING project.

LOCAL COMPUTER
  - It is recommended that you use the code editor ATOM. It is integrated with
    GitHub and contains many of the basic functions needed to change and update
    the code.
    -  atom.io
  Using ATOM to make changes to the code
    - ctrl+shift+p opens up a command window where you can call functions of
      ATOM. You can use Git:Clone, or alt+G+=, to clone a repository to your
      local machine.
    - The features for handling git hub are in the Git and GitHub tabs.
    - In the git tab, it will track the changes that you have made to your
      local version of the code in the Unstaged Changes window.
    - When you want to upload the change you have made:
      - Click the Stage All button in the Unstaged Changes window
      - Write a Commit message with the title of the change and a description
        following a return character.
      - Click Commit to master
      - Finally, click Push at the bottom of the tab.
    - At the bottom of the tab, there is a branch button and an action button.
      - The branch button can access different branches in the project and can
        create new ones.
      - The action button will display the suggested action. Right click on it
        to see all actions.

CAMPUS CLUSTER
  - Download repository to cluster: git clone https://github.com/pcarzon/ICCING.git
  - Update repository: git pull
  - Upload changes to github:
    - git add .
    - git commit -m "Commit Message here"
    - git push
  - See branches of project: git branch -a
  - Switch to different branch: git checkout Test
  - Merge a branch to the main: git merge branch Test
