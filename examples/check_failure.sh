# function to test the exit status of a job
check_failure () {
    # usage: check_failure $?
    if test $1 != 0
    then
        $ECHO "Error condition encountered during test: exit status = $1"
        $ECHO "Aborting"
        exit 1
    fi
}
