#! /usr/bin/python -tt

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from getpass import getpass
from smtplib import SMTP

"""
This script helps to check that the SMTP_HOST (see below) would accept STARTTLS
command, and if LOCAL_HOST is acceptable for it, would check the requested user
name and password would allow to send e-mail through it.
"""



SMTP_HOST = 'smtp.gmail.com:587'
LOCAL_HOST = ''


def main():
    """
    entry point
    """

    server = SMTP(SMTP_HOST)

#    server.ehlo()
#    print(server.ehlo())
 
    server.starttls()

    print(server.ehlo(LOCAL_HOST))

    user = raw_input('user: ')
    password = getpass('password: ')

    print(server.login(user, password))

    fromaddr = 'testfarmqef@gmail.com'
    toaddrs  = 'samuel.pon@gmail.com'
    msg = "\r\n".join([
      "From: testfarmqef@gmail.com",
      "To: samuel.pon@gmail.com",
      "Subject: Buildbot",
      "",
      "Why, oh why"
      ])


    server.sendmail(fromaddr, toaddrs, msg) 
    server.close()

if __name__ == '__main__':
    main()

# vim:ts=4:sw=4:et:tw=80
