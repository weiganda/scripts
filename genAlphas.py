#!/usr/bin/env python

import argparse

def generate_alphas(alpha_min, alpha_max, N):
    r = (alpha_max / alpha_min) ** (1 / (N - 1))
    alphas = [alpha_min * r**i for i in range(N)]
    return alphas

def main():
    parser = argparse.ArgumentParser(description="Generate geometric alphas for basis sets.")
    parser.add_argument('-min', type=float, required=True, help='Minimum (most diffuse) alpha value')
    parser.add_argument('-max', type=float, required=True, help='Maximum (most contracted) alpha value')
    parser.add_argument('-num', type=int, required=True, help='Number of alphas to generate')

    args = parser.parse_args()

    alphas = generate_alphas(args.min, args.max, args.num)
    for a in alphas:
        print(f"{a:.8E}")

if __name__ == "__main__":
    main()
