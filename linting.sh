#! /usr/bin/env bash

isort src tests demo

black src tests demo

flake8 src demo