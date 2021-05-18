# Once submodules are updated, use this script to update submodule_commit_hash_records based on git repo database
git ls-tree HEAD . | grep "^160000" | awk '{print $3, $4}' > submodule_commit_hash_records
