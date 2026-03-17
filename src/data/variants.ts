export interface Variant {
  aa: string;
  score: number;
  misha?: boolean;
}

export const variants: Variant[] = [
  { aa: 'A', score: 0.9016, misha: true },
  { aa: 'C', score: 0.9984 }, { aa: 'D', score: 0.9483 },
  { aa: 'F', score: 0.9992 }, { aa: 'G', score: 0.9191 },
  { aa: 'H', score: 0.9927 }, { aa: 'I', score: 0.9923 },
  { aa: 'K', score: 0.9272 }, { aa: 'L', score: 0.9929 },
  { aa: 'M', score: 0.9909 }, { aa: 'N', score: 0.9822 },
  { aa: 'P', score: 0.9985 }, { aa: 'Q', score: 0.846 },
  { aa: 'R', score: 0.9634 }, { aa: 'S', score: 0.9433 },
  { aa: 'T', score: 0.9666 }, { aa: 'V', score: 0.9664 },
  { aa: 'W', score: 0.9997 }, { aa: 'Y', score: 0.9981 },
];
