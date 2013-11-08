#!/bin/sh
# sync windows installers to repository

# function to selectively remove older packages
# and make a symlink to the latest version
prune_exe () {
  dir=$1
  ref=$(date +%s)

  for bit in 32 64
  do \
    for var in serial mpich2
    do \
      old=999999999
      for exe in ${dir}/qe4win-*svn20[0-9][0-9]*-${bit}bit-${var}.exe
      do \
        [ -f ${exe} ] || continue
        # re-set symbolic link to latest entry
        p=$(echo ${exe} | sed -e 's@^\(.*/\)\(qe4win-.*\)\(20[0-9][0-9]\+\).*\(\.exe\)$@\1@')
        r=$(echo ${exe} | sed -e 's@^\(.*/\)\(qe4win-.*\)\(20[0-9][0-9]\+\).*\(\.exe\)$@\2@')
        t=$(echo ${exe} | sed -e 's@^\(.*/\)\(qe4win-.*\)\(20[0-9][0-9]\+\).*\(\.exe\)$@\3@')
        e=$(echo ${exe} | sed -e 's@^\(.*/\)\(qe4win-.*\)\(20[0-9][0-9]\+\).*\(\.exe\)$@\4@')

        # compute age difference in days
        y=$(echo ${t} | cut -c 1-4)
        m=$(echo ${t} | cut -c 5-6)
        d=$(echo ${t} | cut -c 7-8)
        s=$(date +%s -d "${m}/${d}/${y}")
        age=$(expr \( $ref - $s \) / 86400)

        if [ $age -lt $old ] 
        then
          old=$age
          sym="${r}${t}-${bit}bit-${var}${e}"
          sto="${p}${r}-latest-${bit}bit-${var}${e}"
        fi

        # NOTE: to simplify the math, for the following we
        # define one month to have 28 days and a year of 
        # 12 months to have correspondingly only 336 days.

        # after about one year we keep only one per year
        unset tmp
        if [ $age -gt 336 ]
        then
          y=$(expr $age / 336)
          eval tmp=\$year${bit}${y}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          else
              echo "first in year $y $exe"
          fi
          export year${bit}${y}=1
        fi

        # after about three months we keep only one per month
        unset tmp
        if [ $age -gt 84 ] && [ $age -lt 336 ]
        then
          m=$(expr $age / 28)
          eval tmp=\$month${bit}${m}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          else
              echo "first in month $m $exe"
          fi
          export month${bit}${m}=1
        fi

        # after one week we keep only one per week.
        unset tmp
        if [ $age -gt 7 ] && [ $age -lt 84 ]
        then
          w=$(expr $age / 7)
          eval tmp=\$week${bit}${w}
          if [ -n "$tmp" ]
          then
              rm -vf ${exe}
          else
              echo "first in week $w $exe"
          fi
          export week${bit}${w}=1
        fi
      done
      rm -f ${sto}
      ln -s ${sym} ${sto}
    done
  done
}

if [ -n "${QE4WIN_REPO_USER}" ] \
&& [ -n "${QE4WIN_REPO_HOST}" ] \
&& [ -n "${QE4WIN_REPO_DIR}" ]
then
    pushd ${HOME}/mingw-cross

    prune_exe qe4win

    rsync -arpv --delete qe4win/ \
        ${QE4WIN_REPO_USER}@${QE4WIN_REPO_HOST}:${QE4WIN_REPO_DIR}/

    ssh ${QE4WIN_REPO_USER}@${QE4WIN_REPO_HOST} "cd ${QE4WIN_REPO_DIR}/../; ./mkhtmlindex.sh qe4win"
    popd
else
    cat <<EOF

Required environment variables to determine the target repository
account, server and location are not fully configured.
QE4WIN_REPO_USER=${QE4WIN_REPO_USER}
QE4WIN_REPO_HOST=${QE4WIN_REPO_HOST}
QE4WIN_REPO_DIR=${QE4WIN_REPO_DIR}

EOF
    exit 1
fi
