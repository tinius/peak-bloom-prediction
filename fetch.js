import {JSDOM } from "jsdom"
import { stringify } from 'csv-stringify/sync'
import fs from 'fs'
import _ from 'lodash'

const indexToDate = (monthIndex, dayIndex) =>
  new Date(Date.UTC(2026, monthIndex + 1, dayIndex + 1));

const urls = ['https://www.accuweather.com/en/us/washington/20006/february-weather/327659?year=2026', 'https://www.accuweather.com/en/us/washington/20006/march-weather/327659?year=2026', 'https://www.accuweather.com/en/us/washington/20006/april-weather/327659?year=2026']

const fToC = f => (f - 32)/1.8

const pseq = (arr, lambda) => {

	return arr.reduce( (agg, cur, i, arr) => {

		return agg.then(res => lambda(cur, i, arr).then( res2 => res.concat(res2)))

	}, Promise.resolve([]) )

}

pseq(urls, (url, monthIndex) => {

    return fetch(url, {
  "headers": {
    "accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
    "accept-language": "en-GB,en;q=0.9,de-AT;q=0.8,de;q=0.7,tr-TR;q=0.6,tr;q=0.5,es-CO;q=0.4,es;q=0.3,en-US;q=0.2",
    "priority": "u=0, i",
    "sec-ch-ua": "\"Chromium\";v=\"142\", \"Google Chrome\";v=\"142\", \"Not_A Brand\";v=\"99\"",
    "sec-ch-ua-mobile": "?0",
    "sec-ch-ua-platform": "\"macOS\"",
    "sec-fetch-dest": "document",
    "sec-fetch-mode": "navigate",
    "sec-fetch-site": "none",
    "sec-fetch-user": "?1",
    "upgrade-insecure-requests": "1",
    },
  "body": null,
  "method": "GET"
}).then(resp => resp.text())
.then(str => {
    
    const dom = new JSDOM(str)

    const data = Array.from(dom.window.document.querySelectorAll('.temp'))
        .map((el, i) => {
            return {
                date : indexToDate(monthIndex, i).toISOString().slice(0, 10),
                tmax : Number(el.querySelector('.high').textContent.trim().slice(0, -1)),
                tmin : Number(el.querySelector('.low').textContent.trim().slice(0, -1)),
            }
        })
        .map(row => {
            return {...row, tavg : (row.tmin+row.tmax)/2}
        })

        .map((row, i, arr) => {

            if(arr[0].tavg > 0) {
                // we are looking at F

                return { ...row, tavg : fToC(row.tavg) }

            }

            return row

        })

    return data

})



}).then((allData) => {
    console.log(allData)

    const deduped = _.uniqBy(allData, 'date')
    console.log(deduped)

    fs.writeFileSync(`forecast_temps/accuweather_dc_${(new Date).toISOString().slice(0, 10)}.csv`, stringify(deduped, { header : true }))
})