#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: async_version.py
@time: 2021/1/22 3:03 PM
"""
import asyncio
from functools import wraps
from time import strftime


def decorate(func):
    @wraps(func)
    async def wrapper(*args, **kwargs):
        print(f'{strftime("[%H:%M:%S]")} this is a wrapper')
        res = await func(*args, **kwargs)
        return res

    return wrapper


@decorate
async def tests(n):
    await asyncio.sleep(n / 2)
    print(f'{strftime("[%H:%M:%S]")} this is a test')
    return n


async def main():
    tasks = [asyncio.create_task(tests(i)) for i in range(5, 10)]
    res = await asyncio.gather(*tasks)
    return res


if __name__ == "__main__":
    print(asyncio.run(main()))
