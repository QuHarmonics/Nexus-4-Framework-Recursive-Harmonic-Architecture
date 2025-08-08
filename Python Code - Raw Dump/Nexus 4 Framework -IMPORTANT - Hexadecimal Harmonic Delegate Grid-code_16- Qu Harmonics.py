import React from "react";
import { useEffect, useRef } from "react";

export default function SHA3Spiral() {
  const canvasRef = useRef(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;

    // SHA constants from cube roots of first 64 primes
    const primes = [
      2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
      59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
      127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
      191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
      257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317
    ];

    const K = primes.map(p => Math.floor((2 ** 32) * (p ** (1 / 3) % 1)));
    const ratios = K.slice(1).map((k, i) => k / K[i]);

    ctx.translate(canvas.width / 2, canvas.height / 2);
    ctx.beginPath();

    let angle = 0;
    let radius = 5;

    for (let i = 0; i < ratios.length; i++) {
      angle += ratios[i] * 0.1;
      radius += 2;
      const x = radius * Math.cos(angle);
      const y = radius * Math.sin(angle);

      ctx.fillStyle = `hsl(${(angle * 180) / Math.PI % 360}, 80%, 60%)`;
      ctx.beginPath();
      ctx.arc(x, y, 4, 0, 2 * Math.PI);
      ctx.fill();
    }
  }, []);

  return <canvas ref={canvasRef} className="w-full h-screen" />;
} 
