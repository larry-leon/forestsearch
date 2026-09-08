#!/bin/bash
# Portable hard timeout: tmo.sh SECONDS cmd args...   (exit 124 on timeout)
perl -e 'my $t = shift @ARGV; my $pid = fork(); if ($pid == 0) { exec @ARGV or exit 127 } local $SIG{ALRM} = sub { kill "TERM", $pid; sleep 5; kill "KILL", $pid; exit 124 }; alarm $t; waitpid($pid, 0); exit($? >> 8);' "$@"
