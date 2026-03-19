/**
 * Per-component i18n: each component loads its own translation files.
 * Usage in component:
 *   import { useComponentTranslations } from '../../lib/component-i18n';
 *   import en from './i18n/MyComponent.en.json';
 *   import zh from './i18n/MyComponent.zh.json';
 *   const t = useComponentTranslations(lang, { en, zh });
 */
export function useComponentTranslations(
  lang: string,
  translations: Record<string, Record<string, string>>
) {
  const strings = translations[lang] || translations.en || {};
  const fallback = translations.en || {};
  return function t(key: string): string {
    return strings[key] ?? fallback[key] ?? key;
  };
}
