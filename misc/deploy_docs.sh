#/bin/bash

set -e

## run on wenjie's droplets
build_dir=$(pwd)
target_dir=$HOME/wenjie/wenjie-stat.me/static/intsurv
cd $HOME/wenjie/wenjie-stat.me
git checkout -f
git checkout master
git pull origin master
mkdir -p $target_dir
cp -r $build_dir/docs/* $target_dir
tmp_log=.git_status.log
git status > $tmp_log
if egrep -q "modified:[ ]+docs/" $tmp_log
then
    git add -u static/reda/
    git commit -m "deploy intsurv $CI_COMMIT_SHORT_SHA by gitlab-runner"
    git push origin master
else
    printf "The docs was not updated.\n"
fi
rm $tmp_log
