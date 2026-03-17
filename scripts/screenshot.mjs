#!/usr/bin/env node
/**
 * Headless screenshot utility for QEC Search dashboard.
 *
 * Usage:
 *   node scripts/screenshot.mjs [options]
 *
 * Options:
 *   --url        Page URL             (default: http://localhost:5179)
 *   --out        Output PNG path      (default: /tmp/qec-screenshot.png)
 *   --width      Viewport width       (default: 1920)
 *   --height     Viewport height      (default: 1080)
 *   --wait       Wait ms after load   (default: 2000)
 *   --click      CSS selectors to click, comma-separated
 *   --full       Full-page screenshot (default: false)
 */

import { chromium } from 'playwright'

const IGNORED_PATTERNS = [
  /structuredClone/,
  /Zustand/i,
]

function isIgnored(msg) {
  return IGNORED_PATTERNS.some(p => p.test(msg))
}

function parseArgs() {
  const args = process.argv.slice(2)
  const opts = {
    url: 'http://localhost:5179',
    out: '/tmp/qec-screenshot.png',
    width: 1920,
    height: 1080,
    wait: 2000,
    click: '',
    full: false,
  }

  for (let i = 0; i < args.length; i++) {
    const key = args[i].replace(/^--/, '')
    if (key === 'full') {
      opts.full = true
      continue
    }
    const val = args[++i]
    if (key in opts) {
      opts[key] = ['width', 'height', 'wait'].includes(key)
        ? parseInt(val, 10)
        : val
    }
  }
  return opts
}

async function main() {
  const opts = parseArgs()
  const errors = []

  const browser = await chromium.launch({ headless: true })
  const context = await browser.newContext({
    viewport: { width: opts.width, height: opts.height },
    deviceScaleFactor: 2,
  })
  const page = await context.newPage()

  page.on('pageerror', e => {
    if (!isIgnored(e.message)) {
      errors.push(`PAGE ERROR: ${e.message}`)
      console.error(`PAGE ERROR: ${e.message}`)
    }
  })
  page.on('console', msg => {
    if (msg.type() === 'error') {
      const text = msg.text()
      if (!isIgnored(text)) {
        errors.push(`CONSOLE ERROR: ${text}`)
        console.error(`CONSOLE ERROR: ${text}`)
      }
    }
  })

  try {
    await page.goto(opts.url, { waitUntil: 'networkidle', timeout: 15000 })
  } catch {
    await page.goto(opts.url, { waitUntil: 'load', timeout: 15000 })
  }

  const hasViteOverlay = await page.$('vite-error-overlay').then(el => !!el).catch(() => false)
  if (hasViteOverlay) {
    errors.push('VITE OVERLAY detected')
    console.error('VITE OVERLAY detected')
  }

  await page.waitForTimeout(opts.wait)

  if (opts.click) {
    const selectors = opts.click.split(',').map(s => s.trim()).filter(Boolean)
    for (const sel of selectors) {
      try {
        await page.click(sel, { timeout: 5000, force: true })
        await page.waitForTimeout(1000)
      } catch (e) {
        console.error(`Warning: click "${sel}" failed: ${e.message}`)
      }
    }
  }

  await page.screenshot({
    path: opts.out,
    fullPage: opts.full,
    timeout: 60000,
  })

  await browser.close()
  console.log(opts.out)

  if (errors.length > 0) {
    console.error(`\n${errors.length} error(s) detected`)
    process.exit(2)
  }
}

main().catch(e => {
  console.error(e.message)
  process.exit(1)
})
