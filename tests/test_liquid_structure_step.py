#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `liquid_structure_step` package."""

import pytest  # noqa: F401
import liquid_structure_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = liquid_structure_step.LiquidStructure()
    assert str(type(result)) == "<class 'liquid_structure_step.liquid_structure.LiquidStructure'>"
