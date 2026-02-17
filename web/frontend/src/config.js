/**
 * Application configuration
 * Environment variables are loaded from .env files by Vite
 */

// API base URL - defaults to localhost for development
export const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000/api';

// Application version
export const APP_VERSION = import.meta.env.VITE_APP_VERSION || '1.0.0';
