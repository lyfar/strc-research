import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

const papers = defineCollection({
  loader: glob({ pattern: '*.md', base: './papers' }),
  schema: z.object({
    title: z.string(),
    authors: z.array(z.string()),
    journal: z.string(),
    date: z.coerce.date(),
    pubmed_id: z.string().optional(),
    biorxiv_id: z.string().optional(),
    doi: z.string().optional(),
    relevance_type: z.string().optional(),
    relevance_score: z.string().optional(),
    tags: z.array(z.string()),
    status: z.string().default('unread'),
    date_added: z.coerce.date().optional(),
    publication_status: z.enum(['ahead-of-print', 'published']).optional(),
  }),
});

export const collections = { papers };
