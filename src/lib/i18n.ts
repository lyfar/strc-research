import en from '../locales/en.json';
import zh from '../locales/zh.json';
import fr from '../locales/fr.json';
import es from '../locales/es.json';
import ja from '../locales/ja.json';
import ru from '../locales/ru.json';

const locales: Record<string, Record<string, string>> = { en, zh, fr, es, ja, ru };

export function useTranslations(lang: string) {
  const strings = locales[lang] || locales.en;
  return function t(key: string): string {
    return strings[key] ?? locales.en[key] ?? key;
  };
}

export const languages = ['en', 'zh', 'fr', 'es', 'ja', 'ru'] as const;
export type Lang = typeof languages[number];
